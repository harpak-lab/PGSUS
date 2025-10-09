import pandas as pd 
import numpy as np 
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import seaborn as sns
import os
import seaborn as sns
import statsmodels.api as sm
import warnings
import argparse
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
warnings.filterwarnings('ignore')

class main_figures(object):

    def __init__(self, analyses, label_dict):

        self.analyses = analyses
        self.label_dict = label_dict
        self.label_dict['mean'] = 'Mean Value'
        self.label_dict['median'] = 'Median'
        self.label_dict[0] = 'Rank mean'

    def run(self):
        print(self.analyses)
        if 'alpha_pillars' in self.analyses:
            self.plot_alphapillars()
        if 'decomps' in self.analyses:
            self.plot_decomps()
        if 'insets' in self.analyses:
            self.plot_insets()
        if 'stratification_comparison' in self.analyses:
            self.plot_stratification_comparison()
        if 'decomps_nondirect' in self.analyses:
            self.plot_decomps_nondirect()
    
    def plot_stratification_comparison(self):
        analysis = '1kg.all'
        thresh = '1e-05'
        trait_types ={'alcohol_intake_freq':'behavioral','household_income':'behavioral','neuroticism_score':'behavioral','overall_health':'behavioral','pack_years_smoking':'behavioral','years_schooling':'behavioral', \
        'birth_weight':'anthropometric','bmi':'anthropometric','height':'anthropometric','hip_circ':'anthropometric','waist_circ':'anthropometric', \
        'dbp':'other','fvc':'other','pulse_rate':'other','hand_grip_strength':'other','hair_color':'other','skin_color':'other'}

        wc_nopcs = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
        wc_ukb_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
        wc_1kg_pcs_only = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.' + analysis + '.pcs.only.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
        wc_ukb_and_1kg_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.ukb.and.' + analysis + '.pcs.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
        bolt_nopcs = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
        bolt_wpcs = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')

        wc_nopcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.sad.pc.count.6.mat.txt', sep = '\t').set_index('Unnamed: 0')
        wc_ukb_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.sad.pc.count.6.mat.txt', sep = '\t').set_index('Unnamed: 0')
        wc_1kg_pcs_only_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.' + analysis + '.pcs.only.sps23.v2.sad.pc.count.6.mat.txt',sep = '\t').set_index('Unnamed: 0')
        wc_ukb_and_1kg_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.ukb.and.' + analysis + '.pcs.sps23.v2.sad.pc.count.6.mat.txt',sep = '\t').set_index('Unnamed: 0')
        bolt_nopcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.sad.pc.count.6.mat.txt',sep = '\t').set_index('Unnamed: 0')
        bolt_wpcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.sad.pc.count.6.mat.txt',sep = '\t').set_index('Unnamed: 0')

        def rand_jitter(arr,nanalysis,counter): 
            return arr - (0.9/nanalysis)*counter + 0.5

        methodlabels = ['No Adjustment','GWAS sample PCs','Prediction sample PCs','Both sets of PCs','Linear Mixed Model\n(LMM)','LMM with GWAS\nsample PCs']
        fig,ax = plt.subplots(2,3,figsize = (10,10), gridspec_kw = {'height_ratios' :[2/19.,17/19.],'width_ratios' :[4.5,1,4.5]})
        temp = wc_ukb_pcs.sort_values(by=thresh)
        temp = temp[[thresh]]
        temp.columns = ['wc.ukb']
        temp['wc.nopcs']=wc_nopcs.loc[temp.index.tolist(),thresh]
        temp['wc.1kg.pcs']=wc_1kg_pcs_only.loc[temp.index.tolist(),thresh]
        temp['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs.loc[temp.index.tolist(),thresh]
        temp['bolt.nopcs']=bolt_nopcs.loc[temp.index.tolist(),thresh]
        temp['bolt.wpcs']=bolt_wpcs.loc[temp.index.tolist(),thresh]
        temp['trait_type'] = temp.index.to_series().map(trait_types)                
        temp = temp.drop(columns=['trait_type'])
        temp_median = pd.DataFrame(temp.median(axis=0),columns = ['median'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
        temp_mean = pd.DataFrame(temp.mean(axis=0),columns = ['mean'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
        temp_groups = temp_median.append(temp_mean)
        temp = temp[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        temp = temp.sort_values(by=['wc.ukb'])
        temp = temp.reindex(index=temp.index[::-1])
        ordering =temp.index.tolist()
        temp = temp.loc[ordering]
        temp = temp[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        temp_rank = temp.rank(axis=1, method = 'min')-1
        temp_rank = temp_rank[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        rank_mean = pd.DataFrame(temp_rank.mean(axis = 0)+1).T

        rank_mean = rank_mean.append(temp_mean)
        rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
        rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]

        sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
        sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
        ax[0,2].xaxis.set_label_position('top')
        ax[0,2].xaxis.tick_top()
        ax[0,2].set_xlabel('Isotropic inflation',fontsize = 16)
        ax[0,2].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
        ax[0,2].set_yticklabels(['' for trait in rank_mean.index.tolist()])
        ax[0,2].yaxis.set_tick_params(labelleft=False)
        ax[0,2].set_yticks([])
        ax[0,2].set_ylabel('')
        temp_plot = temp.rank(axis=1, method = 'min')-1
        sns.heatmap(temp_plot,annot=False,ax=ax[1,2],cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
        sns.heatmap(temp_plot,annot=temp, ax = ax[1,2],cmap='Oranges',fmt=".2f", cbar=False)
        ax[1,2].yaxis.set_tick_params(labelleft=False)
        ax[1,2].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
        ax[1,2].set_yticks([])
        ax[1,2].set_ylabel('')
        ax[1,2].set_xticks([])
        alphas = temp
        ax[0,1].axis('off')

        for z,trait in enumerate(rank_mean_rank.index.tolist()[::-1]):
            ax[0,1].text(0.5,(z/2.)+0.25,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
        ax[1,1].axis('off')
        ax[1,1].set_xlim(ax[0,1].get_xlim())
        for z,trait in enumerate(temp.index.tolist()[::-1]):
            ax[1,1].text(0.5,z/17.+0.03,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
        temp = wc_ukb_pcs_sads.sort_values(by=thresh)
        temp = temp[[thresh]]
        temp.columns = ['wc.ukb']
        temp['wc.nopcs']=wc_nopcs_sads.loc[temp.index.tolist(),thresh]
        temp['wc.1kg.pcs']=wc_1kg_pcs_only_sads.loc[temp.index.tolist(),thresh]
        temp['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs_sads.loc[temp.index.tolist(),thresh]
        temp['bolt.nopcs']=bolt_nopcs_sads.loc[temp.index.tolist(),thresh]
        temp['bolt.wpcs']=bolt_wpcs_sads.loc[temp.index.tolist(),thresh]
        temp['trait_type'] = temp.index.to_series().map(trait_types)
        temp = temp.drop(columns=['trait_type'])
        temp_groups = pd.DataFrame(temp.median(axis=0),columns = ['all'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
        temp_groups = temp_groups.reindex(index=temp_groups.index[::-1])
        temp_groups = temp_groups[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        temp = temp.loc[ordering]
        temp = temp[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        temp_rank = temp.rank(axis=1, method = 'min')
        temp_rank-=1
        rank_mean = pd.DataFrame(temp_rank.mean(axis = 0)+1).T
        rank_mean = rank_mean.append(temp_mean)
        rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
        rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
        rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]
        sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
        sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
        ax[0,0].xaxis.set_label_position('top')
        ax[0,0].xaxis.tick_top()
        ax[0,0].set_xlabel('Number of significant PC-wise\nSAD components in top 6 PCs',fontsize = 16)
        ax[0,0].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
        ax[0,0].set_yticklabels(['' for trait in rank_mean.index.tolist()])
        ax[0,0].yaxis.set_tick_params(labelleft=False)
        ax[0,0].set_yticks([])
        ax[0,0].set_ylabel('')
        temp_plot = temp.rank(axis=1, method = 'min')-1
        temp = temp.astype(int)
        sns.heatmap(temp_plot,annot=False,ax=ax[1,0],cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
        sns.heatmap(temp_plot,annot=temp, ax = ax[1,0],cmap='Oranges',fmt="d", cbar=False)
        ax[1,0].yaxis.set_tick_params(labelleft=False)
        ax[1,0].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
        ax[1,0].set_yticks([])
        ax[1,0].set_ylabel('')
        ax[1,0].set_xticks([])
        plt.suptitle('Comparison of SAD variance among methods to adjust for population structure',fontsize=16)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.14)
        sns.despine(bottom=True, top=False, right=True,ax = ax[0,0])
        sns.despine(bottom=True, top=False, right=True,ax = ax[0,2])
        sns.despine(top=True, right=True, bottom=True, ax = ax[1,0])
        sns.despine(top=True, right=True, bottom=True, ax = ax[1,2])
        fig.subplots_adjust(top=0.8)
        plt.savefig('../figures/main_text/strat.control.comparison.' + analysis + '.pval.' + thresh + '.pdf')

    def plot_decomps(self):
        pcs=6
        palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41', 'nondirect':'#8a461b'}
        fig, ax = plt.subplots(nrows = 2, ncols = 4, figsize = ((4*pcs*0.8+1),10))
        
        #top left panel
        df = pd.read_csv('../cache/component_inputs/nondirect/giant/giant_height.1kg.eur.block.permutation.stats.pval.1.0.txt', sep = '\t')
        col_nums = np.array([int(i+1) for i in range(pcs)])
        df = df.set_index('Unnamed: 0')
        ax[0,0].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[0,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,0].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,0].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,0].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[0,0].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[0,0].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[0,0].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[0,0].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[0,0].get_ylim()[1]
        lower_y = ax[0,0].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,0].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,0].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,0].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,0].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,0].set_xlim(0.5, pcs+0.5)
        
        #top middle panel
        df = pd.read_csv('../cache/component_inputs/nondirect/giant2022/giant.eur.1kg.eur.block.permutation.stats.pval.1.0.txt', sep = '\t')
        col_nums = np.array([int(i+1) for i in range(pcs)])
        df = df.set_index('Unnamed: 0')
        ax[0,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[0,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,1].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,1].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,1].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[0,1].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[0,1].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[0,1].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[0,1].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[0,1].get_ylim()[1]
        lower_y = ax[0,1].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,1].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,1].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,1].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,1].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,1].set_xlim(0.5, pcs+0.5)

        #okbay
        df = pd.read_csv('../cache/component_inputs/nondirect/okbay2022/okbay.1kg.eur.block.permutation.stats.pval.0.001.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[0,2].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[0,2].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,2].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,2].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,2].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[0,2].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[0,2].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[0,2].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[0,2].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[0,2].get_ylim()[1]
        lower_y = ax[0,2].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,2].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,2].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,2].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,2].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,2].set_xlim(0.5, pcs+0.5)


        # bmi plink.wc ukb - 1kg all
        df = pd.read_csv('../cache/component_inputs/nondirect/wc/plink.wc.1kg.all.sps23.bmi.aperm.1K.to.1M.block.permutation.stats.pval.1e-05.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[1,0].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[1,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,0].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,0].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,0].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[1,0].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[1,0].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[1,0].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[1,0].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[1,0].get_ylim()[1]
        lower_y = ax[1,0].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,0].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,0].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,0].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,0].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,0].set_xlim(0.5, pcs+0.5)


        #top right forced vital capacity plink.wc ukb
        df = pd.read_csv('../cache/component_inputs/nondirect/wc/plink.wc.1kg.eur.sps23.fvc.aperm.1K.to.1M.block.permutation.stats.pval.0.001.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])
        
        ax[1,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[1,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,1].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,1].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,1].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')

        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[1,1].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[1,1].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[1,1].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[1,1].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[1,1].get_ylim()[1]
        lower_y = ax[1,1].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,1].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,1].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,1].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,1].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,1].set_xlim(0.5, pcs+0.5)


        #EA plink.wc ukb
        df = pd.read_csv('../cache/component_inputs/nondirect/wc/plink.wc.1kg.eur.sps23.years_schooling.aperm.1K.to.1M.block.permutation.stats.pval.1.0.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[1,2].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[1,2].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,2].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,2].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,2].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[1,2].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[1,2].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[1,2].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[1,2].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[1,2].get_ylim()[1]
        lower_y = ax[1,2].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,2].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,2].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,2].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,2].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,2].set_xlim(0.5, pcs+0.5)

        #bottom center
        df = pd.read_csv('../cache/component_inputs/nondirect/akbari2024/akbari2024.overall_health.block.permutation.stats.pval.1.0.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])
        
        ax[1,3].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,3].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,3].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,3].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

        ax[1,3].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,3].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,3].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,3].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
        df = df[df.columns[:6]]
        df.columns = [i for i in range(6)]
        direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
        sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
        covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
        nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

        if len(direct_sig_cols) != 0:
            ax[1,3].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(sad_sig_cols) != 0:
            ax[1,3].scatter(col_nums[sad_sig_cols]-0.125, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(covar_sig_cols) != 0:
            ax[1,3].scatter(col_nums[covar_sig_cols]+0.125, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
        if len(nondirect_sig_cols) != 0:
            ax[1,3].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

        upper_y = ax[1,3].get_ylim()[1]
        lower_y = ax[1,3].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,3].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,3].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,3].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,3].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,3].set_xlim(0.5, pcs+0.5)

        ax[0,3].axis('off')
        sns.despine()
        plt.tight_layout()
        plt.savefig('../figures/main_text/fig.significant.component.main.examples.pdf')
        
    def plot_insets(self):
        #make the pc space for each of the insets
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        superpop_color_dict = {'AFR':'#e7298a','EUR':'#9e771b','EAS':'#7570b3','AMR':'#d95f02','SAS':'#027cd9'}
        eur_pop_color_dict = {'TSI':'#027cd9','FIN':'#9e771b','CEU':'#7570b3','GBR':'#d95f02','IBS':'#e7298a'}
        #giant - panel A
        labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
        indiv_pop = dict(zip(labeldata['sample'],labeldata['color']))
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/giant.height.eur.flashpcs',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        # eigenvals = pd.read_csv('../cache/inset_data/giant.height.eur.eigenvals',sep = '\t',header = None)
        # eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (%)' )
        ax.set_ylabel('PC2 (%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/giant.height.eur.PC1.v.PC2.pdf')
        plt.clf()
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_pop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/years_schooling.1kg.eur.pc',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        eigenvals = pd.read_csv('../cache/1kg_pc_data/years_schooling.1kg.eur.eigenvalues',sep = '\t',header = None, dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/years_schooling.eur.PC1.v.PC2.pdf')
        plt.clf()
        #ukb fvc - panel b
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_pop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/fvc.1kg.eur.pc',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        eigenvals = pd.read_csv('../cache/1kg_pc_data/fvc.1kg.eur.eigenvalues',sep = '\t',header = None, dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/fvc.eur.PC1.v.PC2.pdf')
        plt.clf()
        #bolt-lmm household income - panel C
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['super_pop'].map(superpop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/bmi.1kg.all.pc',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        eigenvals = pd.read_csv('../cache/1kg_pc_data/bmi.1kg.all.eigenvalues',sep = '\t',header = None, dtype=float)

        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/bmi.all.PC1.v.PC2.pdf')
        plt.clf()
        #ukb waist circ - panel D
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['super_pop'].map(superpop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/waist_circ.1kg.all.pc',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        eigenvals = pd.read_csv('../cache/1kg_pc_data/waist_circ.1kg.all.eigenvalues',sep = '\t',header = None,dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/waist_circ.all.PC1.v.PC2.pdf')
        plt.clf()
        #plot okbay
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/okbay2022.pcs',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        eigenvals = pd.read_csv('../cache/1kg_pc_data/okbay2022.values',sep = '\t',header = None,dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/okbay.eur.PC1.v.PC2.pdf')
        plt.clf()
        #plot yengo
        labeldata = pd.read_csv('../cache/1kg_pc_data/1kg_poplabel_map.txt', sep = '\t')
        labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})
        indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
        gbr = labeldata[labeldata['pop']=='GBR']
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/yengo2022.pcs',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
        gbr = df[df['IID'].isin(gbr['IID'].tolist())]
        eigenvals = pd.read_csv('../cache/1kg_pc_data/yengo2022.values',sep = '\t',header = None,dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/yengo.eur.PC1.v.PC2.pdf')
        plt.clf()

        # plot akbari
        ancient_pops_color_dict = {'AN':'#4e79a7','BA':'#f28e2b','EN':'#e15759','H':'#76b7b2','M':'#59a14f','S':'#edc948'}
        labeldata = pd.read_csv('../cache/1kg_pc_data/le2022-sample-metadata.txt', sep = '\t')

        labeldata['color'] = labeldata['Epoch'].map(ancient_pops_color_dict)
        labeldata = labeldata.rename(columns = {'sample':'IID'})

        indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
        df = pd.read_csv('../cache/1kg_pc_data/le-aadr-pca.v2.pcs',sep = '\t')
        df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')

        eigenvals = pd.read_csv('../cache/1kg_pc_data/le-aadr-pca.v2.eigenvalues',sep = '\t',header = None,dtype=float)
        eigenvals = (eigenvals/np.sum(eigenvals))*100
        ax.scatter(df['PC1'],df['PC2'], s=7, color = df['color'])
        ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
        ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[2][0],2)) + '%)')
        plt.tight_layout()
        plt.savefig('../figures/main_text/akbari2024.overall_health.PC1.v.PC2.pdf')
        plt.clf()

    def plot_alphapillars(self):
        label = '1kg.all'
        alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
        alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
        
        alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
        alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]

        alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
        alpha_dfplot = alpha_dfplot[['1e-05','1e-08','ycoordinate']]
        print(alpha_dfplot.mean())
        fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (9,5))
        for i,j in enumerate(alpha_dfplot.columns[:2]):                 
            ax[i].vlines(1, ymin = 0, ymax =17.5, color = 'black')
            ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
            ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
            ax[i].title.set_text('p < ' + str(j))
            ax[i].set_axisbelow(True)
            ax[i].set_xlabel('Isotropic inflation factor')
            ax[i].grid(True, linewidth = 1)

        ax[0].set_yticks(alpha_dfplot['ycoordinate'])
        ax[0].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
        
        ax[1].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
        ax[1].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
        ax[0].set_xticks([x for x in range(1,10)])
        ax[0].set_xticklabels([str(x) for x in ['1','2','3','4','5','6','7','8','9']])
        ax[1].set_xticks([x for x in range(1,10)])
        ax[1].set_xticklabels([str(x) for x in ['1','2','3','4','5','6','7','8','9']])

        ax[0].title.set_text(r'Marginal GWAS $p$-value < ' + r'$10^{-5}$')
        ax[1].title.set_text(r'Marginal GWAS $p$-value < ' + r'$10^{-8}$')

        ax[0].set_ylim([0.5,17.5])
        ax[1].set_ylim([0.5,17.5])
        ax[0].set_xlim([1,9.5])
        ax[1].set_xlim([1,9.5])
        sns.despine()
        plt.tight_layout()
        plt.savefig('../figures/main_text/fig.wc.' + label + '.alpha.pillar.pdf')
        plt.clf()

    def plot_decomps_nondirect(self):
        pcs=6
        palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41'}
        fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = ((3*pcs*0.8+1),10))
        df = pd.read_csv('../cache/component_inputs/nondirect_test/giant.1kg.eur.nondirect.test.block.permutation.stats.pval.1.0.txt', sep = '\t')
        col_nums = np.array([int(i+1) for i in range(pcs)])
        df = df.set_index('Unnamed: 0')
        ax[0,0].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,0].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')


        ax[0,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,0].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,0].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,0].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')

        upper_y = ax[0,0].get_ylim()[1]
        lower_y = ax[0,0].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,0].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,0].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,0].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,0].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,0].set_xlim(0.5, pcs+0.5)
        ax[0,0].set_ylabel('Component / Total PGS Variance')
        
        #top middle forced vital capacity plink.wc ukb
        df = pd.read_csv('../cache/component_inputs/nondirect_test/plink.wc.1kg.eur.sps23.fvc.nondirect.test.block.permutation.stats.pval.0.001.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[0,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,1].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')

        ax[0,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,1].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,1].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,1].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')
        
        upper_y = ax[0,1].get_ylim()[1]
        lower_y = ax[0,1].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,1].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,1].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,1].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,1].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,1].set_xlim(0.5, pcs+0.5)
        ax[0,1].set_ylabel('Component / Total PGS Variance')

        #top right bmi plink.wc ukb - 1kg all
        df = pd.read_csv('../cache/component_inputs/nondirect_test/plink.wc.1kg.all.sps23.bmi.aperm.nondirect.test.block.permutation.stats.pval.1e-05.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[0,2].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[0,2].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')

        ax[0,2].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[0,2].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[0,2].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[0,2].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')
        
        upper_y = ax[0,2].get_ylim()[1]
        lower_y = ax[0,2].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[0,2].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[0,2].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[0,2].set_xticks([int(i+1) for i in range(pcs)])
        ax[0,2].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[0,2].set_xlim(0.5, pcs+0.5)
        ax[0,2].set_ylabel('Component / Total PGS Variance')

        #bottom left EA plink.wc ukb
        df = pd.read_csv('../cache/component_inputs/nondirect_test/plink.wc.1kg.eur.sps23.years_schooling.nondirect.test.block.permutation.stats.pval.1.0.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[1,0].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,0].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')

        ax[1,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,0].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,0].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,0].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')
        
        upper_y = ax[1,0].get_ylim()[1]
        lower_y = ax[1,0].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,0].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,0].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,0].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,0].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,0].set_xlim(0.5, pcs+0.5)
        ax[1,0].set_ylabel('Component / Total PGS Variance')


        #bottom middle EA plink.wc ukb
        df = pd.read_csv('../cache/component_inputs/nondirect_test/bolt.wpcs.1kg.eur.years_schooling.nondirect.test.block.permutation.stats.pval.1.0.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[1,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,1].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')

        ax[1,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,1].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,1].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,1].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')
        
        upper_y = ax[1,1].get_ylim()[1]
        lower_y = ax[1,1].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,1].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,1].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,1].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,1].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,1].set_xlim(0.5, pcs+0.5)
        ax[1,1].set_ylabel('Component / Total PGS Variance')

        df = pd.read_csv('../cache/component_inputs/nondirect_test/plink.wc.1kg.all.sps23waist_circ.nondirect.test.block.permutation.stats.pval.1e-08.txt', sep = '\t').set_index('Unnamed: 0')
        col_nums = np.array([int(i+1) for i in range(pcs)])

        ax[1,2].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums-0.125, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums+0.125, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
        ax[1,2].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor='#808080', alpha=0.2, edgecolor = 'none')

        ax[1,2].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
        ax[1,2].scatter(col_nums-0.125, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
        ax[1,2].scatter(col_nums+0.125, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
        ax[1,2].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = '#808080', label='nondirect variance')
        
        upper_y = ax[1,2].get_ylim()[1]
        lower_y = ax[1,2].get_ylim()[0]
        asterisk_y = np.array([upper_y for i in range(df.shape[1])])
        ax[1,2].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
        for i in col_nums:
            ax[1,2].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
        ax[1,2].set_xticks([int(i+1) for i in range(pcs)])
        ax[1,2].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
        ax[1,2].set_xlim(0.5, pcs+0.5)
        ax[1,2].set_ylabel('Component / Total PGS Variance')
        
        sns.despine()
        plt.tight_layout()
        plt.savefig('../figures/main_text/fig.significant.component.main.examples.w.nondirect.pdf')


