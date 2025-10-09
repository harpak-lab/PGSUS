import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 600

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--n-pcs", type=int, dest="n_pcs", required=False, default=6)
parser.add_argument("--output", type=str, default = 'sad_decomposition')
args = parser.parse_args()

pfile = args.input
pcs = args.n_pcs
outstem = args.output
outfile = f"{outstem}.png"

palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41', 'nondirect':'#8a461b'}

# df = pd.read_csv('../cache/component_inputs/giant/' + trait + '/giant.' + label + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t')
df = pd.read_table(pfile, sep="\t", header=0, index_col=0)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (pcs*0.67+1,5))
col_nums = np.array([int(i+1) for i in range(pcs)])
ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:pcs], df.loc['lower0_perm_direct'].to_numpy()[:pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:pcs] > df.loc['lower0_perm_direct'].to_numpy()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].to_numpy()[:pcs], df.loc['lower0_perm_sad'].to_numpy()[:pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:pcs] > df.loc['lower0_perm_sad'].to_numpy()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].to_numpy()[:pcs], df.loc['lower025_perm_covar'].to_numpy()[:pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:pcs] > df.loc['lower025_perm_covar'].to_numpy()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].to_numpy()[:pcs], df.loc['lower025_perm_nondirect'].to_numpy()[:pcs], where=(df.loc['upper975_perm_nondirect'].to_numpy()[:pcs] > df.loc['lower025_perm_nondirect'].to_numpy()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:pcs], color = palette['direct'], label='direct variance')
ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].to_numpy()[:pcs], color = palette['sad'], label='SAD variance')
ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].to_numpy()[:pcs], color = palette['covar'], label='direct-SAD covariance')
ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].to_numpy()[:pcs], color = palette['nondirect'], label='nondirect variance')

df = df[df.columns[:6]]
df.columns = [i for i in range(6)]
direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:pcs] < 0.05] # Two-tailed test
sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:pcs] < 0.05]
covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:pcs] < 0.025] # One-tailed test
nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:pcs] < 0.025]

if len(direct_sig_cols) != 0:
    ax.scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
if len(sad_sig_cols) != 0:
    ax.scatter(col_nums[sad_sig_cols]-0.083, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
if len(covar_sig_cols) != 0:
    ax.scatter(col_nums[covar_sig_cols]+0.083, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
if len(nondirect_sig_cols) != 0:
    ax.scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

upper_y = ax.get_ylim()[1]
lower_y = ax.get_ylim()[0]
asterisk_y = np.array([upper_y for i in range(df.shape[1])])

ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
for i in col_nums:
    ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)

ax.set_xticks([int(i+1) for i in range(pcs)])
ax.set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
ax.set_xlim(0.5, pcs+0.5)
ax.set_ylabel('Component / Total PGS variance')
sns.despine()

plt.tight_layout()
fig = ax.get_figure()
fig.savefig(outfile)

