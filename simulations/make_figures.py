import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats
import os
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--figs',type=str,help='figures to make')
args = parser.parse_args()
figs = [int(x) for x in args.figs.split(",")]

#need to go back on TACC and generate all of the 2+ gwas pc simulations with an ascp of 1

####################################################################
#Generate matrix that is ascp thresholds by number of gwas pcs
#and keep the number of pcs in the gamma estimation constant at 100
#if it already exist, skip it

if 1 in figs:
	if os.path.isfile('cache/weighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt'):
		outdf_weighted = pd.read_csv('cache/weighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt', sep = '\t')
		sedf_weighted = pd.read_csv('cache/weighted.ses.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt', sep = '\t')
	else:
		cutoffs = ['1.0','0.1','0.01','0.001','0.0001','1e-05','1e-06','1e-07','1e-08']
		gwaspcs = [i for i in range(1,21)]
		outdf_weighted = pd.DataFrame(np.zeros((len(gwaspcs), len(cutoffs))), index = gwaspcs, columns = cutoffs)
		sedf_weighted = pd.DataFrame(np.zeros((len(gwaspcs), len(cutoffs))), index = gwaspcs, columns = cutoffs)

		for i in cutoffs:
			for j in gwaspcs:
				temp = np.array(pd.read_csv('gamma_outputs/weighted.gammas.' + str(j) + '.gwaspcs.100.pc.lower.bound.' + i + '.asc.p.txt', sep = '\t', header = None))
				outdf_weighted.loc[j,i] = np.mean(temp)
				sedf_weighted.loc[j,i] = scipy.stats.sem(temp)

		outdf_weighted.to_csv('cache/weighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt',sep = '\t')
		sedf_weighted.to_csv('cache/weighted.ses.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt',sep = '\t')

	####################################################################
	#Repeat for the unweighted gammas estimates, if already made skip it

	if os.path.isfile('cache/unweighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt'):
		outdf_unweighted = pd.read_csv('cache/unweighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt', sep = '\t')
		sedf_unweighted = pd.read_csv('cache/unweighted.ses.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt', sep = '\t')
	else:
		cutoffs = ['1.0','0.1','0.01','0.001','0.0001','1e-05','1e-06','1e-07','1e-08']
		gwaspcs = [i for i in range(1,21)]
		outdf_unweighted = pd.DataFrame(np.zeros((len(gwaspcs), len(cutoffs))), index = gwaspcs, columns = cutoffs)
		sedf_unweighted = pd.DataFrame(np.zeros((len(gwaspcs), len(cutoffs))), index = gwaspcs, columns = cutoffs)

		for i in cutoffs:
			for j in gwaspcs:
				temp = np.array(pd.read_csv('gamma_outputs/unweighted.gammas.' + str(j) + '.gwaspcs.100.pc.lower.bound.' + i + '.asc.p.txt', sep = '\t', header = None))
				outdf_unweighted.loc[j,i] = np.mean(temp)
				sedf_unweighted.loc[j,i] = scipy.stats.sem(temp)

		outdf_unweighted.to_csv('cache/unweighted.gammas.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt',sep = '\t')
		sedf_unweighted.to_csv('cache/unweighted.ses.100.pc.lb.ascp.range.e1_8.gwas.pc.range.1_20.txt',sep = '\t')


	fig, ax = plt.subplots(1,2, figsize = (8,4))
	thresh_dict = {'1.0':'o','0.1':'s','0.001':'*','1e-05':'^','1e-08':'8'}
	color_dict = {'1.0':'#762a83','0.1':'#c2a5cf','0.001':'#d9f0d3','1e-05':'#a6dba0','1e-08':'#1b7837'}

	subset_thresh = ['1.0','0.1','0.001','1e-05','1e-08']
	outdf_weighted_plot = outdf_weighted[subset_thresh]
	outdf_weighted_plot = outdf_weighted_plot.reset_index()
	outdf_weighted_plot['index'] += 1

	ax[0].hlines(1, xmin = 0, xmax = 21, color = 'black', linestyle = '--')
	ax[1].hlines(1, xmin = 0, xmax = 21, color = 'black', linestyle = '--', label = r'$\gamma$')

	for thresh in subset_thresh:
		ax[0].plot(np.array(outdf_weighted_plot['index']),np.array(outdf_weighted_plot[thresh]), marker = thresh_dict[thresh], markersize = 4, linestyle = ':', color = color_dict[thresh])

	outdf_unweighted_plot = outdf_unweighted[subset_thresh]
	outdf_unweighted_plot = outdf_unweighted_plot.reset_index()
	outdf_unweighted_plot['index'] += 1

	for thresh in subset_thresh:
		ax[1].plot(np.array(outdf_unweighted_plot['index']),np.array(outdf_unweighted_plot[thresh]), marker = thresh_dict[thresh], markersize = 4, linestyle = ':', color = color_dict[thresh], label = thresh)


	ax[0].set_title('Weighted estimates')
	ax[1].set_title('Unweighted estimates')

	ax[0].set_ylabel(r'$\frac{1}{\hat{\gamma}}$')
	ax[0].set_xlabel('Number of PCs included as \ncovariate in population GWAS')
	ax[1].set_xlabel('Number of PCs included as \ncovariate in population GWAS')

	ax[1].legend(bbox_to_anchor = [1.4,1], ncol = 1, loc = 'upper right')
	plt.tight_layout()
	sns.despine()
	plt.savefig('figures/vary.gwaspcs.ascp.100lb.pdf')


if 2 in figs:

	if os.path.isfile('cache/weighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt'):
		outdf_weighted = pd.read_csv('cache/weighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt', sep = '\t')
		sedf_weighted = pd.read_csv('cache/weighted.ses.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt', sep = '\t')

	else:
		cutoffs = ['1.0','0.1','0.01','0.001','0.0001','1e-05','1e-06','1e-07','1e-08']
		lbpcs = [i for i in range(10,210,10)]
		outdf_weighted = pd.DataFrame(np.zeros((len(lbpcs), len(cutoffs))), index = lbpcs, columns = cutoffs)
		sedf_weighted = pd.DataFrame(np.zeros((len(lbpcs), len(cutoffs))), index = lbpcs, columns = cutoffs)

		for i in cutoffs:
			for j in lbpcs:
				temp = np.array(pd.read_csv('gamma_outputs/weighted.gammas.5.gwaspcs.' + str(j) + '.pc.lower.bound.' + i + '.asc.p.txt', sep = '\t', header = None))
				outdf_weighted.loc[j,i] = np.mean(temp)
				sedf_weighted.loc[j,i] = scipy.stats.sem(temp)

		outdf_weighted.to_csv('cache/weighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt',sep = '\t')
		sedf_weighted.to_csv('cache/weighted.ses.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt',sep = '\t')

	#repeat for the unweighted
	if os.path.isfile('cache/unweighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt'):
		outdf_unweighted = pd.read_csv('cache/unweighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt', sep = '\t')
		sedf_unweighted = pd.read_csv('cache/unweighted.ses.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt', sep = '\t')

	else:
		cutoffs = ['1.0','0.1','0.01','0.001','0.0001','1e-05','1e-06','1e-07','1e-08']
		lbpcs = [i for i in range(10,210,10)]
		outdf_unweighted = pd.DataFrame(np.zeros((len(lbpcs), len(cutoffs))), index = lbpcs, columns = cutoffs)
		sedf_unweighted = pd.DataFrame(np.zeros((len(lbpcs), len(cutoffs))), index = lbpcs, columns = cutoffs)

		for i in cutoffs:
			for j in lbpcs:
				temp = np.array(pd.read_csv('gamma_outputs/unweighted.gammas.5.gwaspcs.' + str(j) + '.pc.lower.bound.' + i + '.asc.p.txt', sep = '\t', header = None))
				outdf_unweighted.loc[j,i] = np.mean(temp)
				sedf_unweighted.loc[j,i] = scipy.stats.sem(temp)

		outdf_unweighted.to_csv('cache/unweighted.gammas.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt',sep = '\t')
		sedf_unweighted.to_csv('cache/unweighted.ses.pc.lb.range.10_200.ascp.range.e1_8.gwas.pc.5.txt',sep = '\t')

	fig, ax = plt.subplots(1,2, figsize = (8,4))
	thresh_dict = {'1.0':'o','0.1':'s','0.001':'*','1e-05':'^','1e-08':'8'}
	color_dict = {'1.0':'#762a83','0.1':'#c2a5cf','0.001':'#d9f0d3','1e-05':'#a6dba0','1e-08':'#1b7837'}

	# print(outdf_weighted)
	subset_thresh = ['Unnamed: 0','1.0','0.1','0.001','1e-05','1e-08']
	outdf_weighted_plot = outdf_weighted[subset_thresh]
	outdf_weighted_plot = outdf_weighted_plot.reset_index()
	outdf_weighted_plot['index'] += 1

	ax[0].hlines(1, xmin = 0, xmax = 21, color = 'black', linestyle = '--')
	ax[1].hlines(1, xmin = 0, xmax = 21, color = 'black', linestyle = '--', label = r'$\gamma$')

	for thresh in subset_thresh[1:]:
		ax[0].plot(np.array(outdf_weighted_plot.index.tolist()),np.array(outdf_weighted_plot[thresh]), marker = thresh_dict[thresh], markersize = 4, linestyle = ':', color = color_dict[thresh])

	outdf_unweighted_plot = outdf_unweighted[subset_thresh]
	outdf_unweighted_plot = outdf_unweighted_plot.reset_index()
	outdf_unweighted_plot['index'] += 1

	for thresh in subset_thresh[1:]:
		ax[1].plot(np.array(outdf_unweighted_plot.index.tolist()),np.array(outdf_unweighted_plot[thresh]), marker = thresh_dict[thresh], markersize = 4, linestyle = ':', color = color_dict[thresh], label = thresh)

	ax[0].set_xticks(outdf_weighted_plot.index.tolist(), fontsize = 12)
	ax[0].set_xticklabels(outdf_weighted['Unnamed: 0'].tolist(), rotation = 60, fontsize = 8)

	ax[1].set_xticks(outdf_weighted_plot.index.tolist())
	ax[1].set_xticklabels(outdf_weighted['Unnamed: 0'].tolist(), rotation = 60, fontsize = 8)

	ax[0].set_title('Weighted estimates')
	ax[1].set_title('Unweighted estimates')

	ax[0].set_ylabel(r'$\frac{1}{\hat{\gamma}}$')
	ax[0].set_xlabel('Lower bound of PCs included in \nSAD estimation of ' + r'$\hat{\gamma}$')
	ax[1].set_xlabel('Lower bound of PCs included in \nSAD estimation of ' + r'$\hat{\gamma}$')

	ax[1].legend(bbox_to_anchor = [1.4,1], ncol = 1, loc = 'upper right')
	plt.tight_layout()
	sns.despine()
	plt.savefig('figures/vary.sadpcs.ascp.5gwaspcs.pdf')

