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
from scipy.stats import pearsonr
import matplotlib as mpl
warnings.filterwarnings('ignore')

class main_figures(object):

	def __init__(self, analyses, label_dict):

		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		if 'alphapillar' in analyses:
			self.plot_alphapillars()
		if 'decomps' in analyses:
			self.plot_decomps()
		if 'stratification_comparison' in analyses:
			self.plot_stratification_comparison()


	def plot_stratification_comparison(self):
		analysis = '1kg.all'
		trait_types ={'alcohol_intake_freq':'behavioral','household_income':'behavioral','neuroticism_score':'behavioral','overall_health':'behavioral','pack_years_smoking':'behavioral','years_schooling':'behavioral', \
		'birth_weight':'anthropometric','bmi':'anthropometric','height':'anthropometric','hip_circ':'anthropometric','waist_circ':'anthropometric', \
		'dbp':'other','fvc':'other','pulse_rate':'other','hand_grip_strength':'other','hair_color':'other','skin_color':'other'}

		wc_nopcs = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
		wc_ukb_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
		wc_1kg_pcs_only = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.' + analysis + '.pcs.only.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		wc_ukb_and_1kg_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.ukb.and.' + analysis + '.pcs.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		bolt_nopcs = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		bolt_wpcs = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')

		wc_nopcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt', sep = '\t').set_index('Unnamed: 0')
		wc_ukb_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.sad.pc.count.mat.txt', sep = '\t').set_index('Unnamed: 0')
		wc_1kg_pcs_only_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.' + analysis + '.pcs.only.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
		wc_ukb_and_1kg_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.ukb.and.' + analysis + '.pcs.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
		bolt_nopcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
		bolt_wpcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')

		def rand_jitter(arr,nanalysis,counter):	
			return arr - (0.9/nanalysis)*counter + 0.5

		methodlabels = ['None','GWAS sample PCs','Target sample PCs','Both sets of PCs','Linear Mixed Model\n(LMM)','LMM with GWAS\nsample PCs']
		thresh = '1e-05'
		fig,ax = plt.subplots(2,3,figsize = (10,10), gridspec_kw = {'height_ratios' :[1/18.,17/18.],'width_ratios' :[4.5,1,4.5]})
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
		temp_groups = pd.DataFrame(temp.median(axis=0),columns = ['all'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
		temp = temp[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
		temp = temp.sort_values(by=['wc.ukb'])
		temp = temp.reindex(index=temp.index[::-1])
		ordering =temp.index.tolist()
		temp_groups = temp_groups[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
		sns.heatmap(temp_groups,annot=True,ax=ax[0,2],cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
		ax[0,2].xaxis.set_label_position('top')
		ax[0,2].xaxis.tick_top()
		ax[0,2].set_xlabel('Isotropic inflation',fontsize = 16)
		ax[0,2].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
		ax[0,2].set_yticklabels(['' for trait in temp_groups.index.tolist()])
		ax[0,2].yaxis.set_tick_params(labelleft=False)
		ax[0,2].set_yticks([])
		ax[0,2].set_ylabel('')
		sns.heatmap(temp,annot=True,ax=ax[1,2],fmt = '.2f',cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
		ax[1,2].yaxis.set_tick_params(labelleft=False)
		ax[1,2].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
		ax[1,2].set_yticks([])
		ax[1,2].set_ylabel('')
		ax[1,2].set_xticks([])
		alphas = temp
		ax[0,1].axis('off')
		for z,trait in enumerate(temp_groups.index.tolist()[::-1]):
			ax[0,1].text(0.5,z + 0.5,label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
		ax[1,1].axis('off')
		ax[1,1].set_xlim(ax[0,1].get_xlim())
		for z,trait in enumerate(temp.index.tolist()[::-1]):
			ax[1,1].text(0.5,z/17.+0.03,label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
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
		sns.heatmap(temp_groups,annot=True,ax=ax[0,0],cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
		ax[0,0].xaxis.set_label_position('top')
		ax[0,0].xaxis.tick_top()
		ax[0,0].set_xlabel('Number of significant PC-wise\nSAD components in first 20 PCs',fontsize = 16)
		ax[0,0].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
		ax[0,0].set_yticklabels(['' for trait in temp_groups.index.tolist()])
		ax[0,0].yaxis.set_tick_params(labelleft=False)
		ax[0,0].set_yticks([])
		ax[0,0].set_ylabel('')
		sns.heatmap(temp,annot=True,ax=ax[1,0],cmap='Oranges',vmin = temp.values.min(),vmax = temp.values.max(),cbar = False)
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
		plt.savefig('../figures/stratification_control_comp/strat.control.comparison.' + analysis + '.pval.' + thresh + '.pdf')

	def plot_decomps(self):
		pcs=6
		palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41'}
		fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = ((2*pcs*0.8+1),10))
		df = pd.read_csv('../cache/component_inputs/giant/giant_height_rescaled/giant.1kg.eur.block.permutation.stats.pval.1.0.txt', sep = '\t')
		col_nums = np.array([int(i+1) for i in range(pcs)])
		df = df.set_index('Unnamed: 0')
		ax[0,0].fill_between(col_nums-0.25, df.loc['upper975_perm_direct'].tolist()[:pcs], df.loc['lower025_perm_direct'].tolist()[:pcs], where=(df.loc['upper975_perm_direct'].tolist()[:pcs] > df.loc['lower025_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
		ax[0,0].fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
		ax[0,0].fill_between(col_nums, df.loc['upper975_perm_sad'].tolist()[:pcs], df.loc['lower025_perm_sad'].tolist()[:pcs], where=(df.loc['upper975_perm_sad'].tolist()[:pcs] > df.loc['lower025_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
		ax[0,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
		ax[0,0].scatter(col_nums, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
		ax[0,0].scatter(col_nums+0.25, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD variance')
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
		#top right forced vital capacity plink.wc ukb
		df = pd.read_csv('../cache/component_inputs/wc/plink.wc.1kg.eur.sps23.fvc.aperm.1K.to.1M.block.permutation.stats.pval.0.001.txt', sep = '\t').set_index('Unnamed: 0')
		col_nums = np.array([int(i+1) for i in range(pcs)])
		ax[0,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
		ax[0,1].fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
		ax[0,1].fill_between(col_nums, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
		ax[0,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
		ax[0,1].scatter(col_nums, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
		ax[0,1].scatter(col_nums+0.25, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD variance')
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
		#bottom left household income bolt-lmm
		df = pd.read_csv('../cache/component_inputs/bolt/bolt.nopcs.1kg.eur.household_income.block.permutation.stats.pval.1.0.txt', sep = '\t').set_index('Unnamed: 0')
		col_nums = np.array([int(i+1) for i in range(pcs)])
		ax[1,0].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
		ax[1,0].fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
		ax[1,0].fill_between(col_nums, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
		ax[1,0].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
		ax[1,0].scatter(col_nums, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
		ax[1,0].scatter(col_nums+0.25, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD variance')
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
		#bottom right waist circ ukb plink.wc
		df = pd.read_csv('../cache/component_inputs/wc/plink.wc.1kg.all.sps23.waist_circ.aperm.1K.to.1M.block.permutation.stats.pval.1e-08.txt', sep = '\t').set_index('Unnamed: 0')
		col_nums = np.array([int(i+1) for i in range(pcs)])
		ax[1,1].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
		ax[1,1].fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
		ax[1,1].fill_between(col_nums, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
		ax[1,1].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
		ax[1,1].scatter(col_nums, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
		ax[1,1].scatter(col_nums+0.25, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD variance')
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
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/main_text/fig.significant.component.main.examples.pdf')
		#make the pc space for each of the insets
		labeldata = pd.read_csv('../cache/inset_data/1kg_poplabel_map.txt', sep = '\t')
		superpop_color_dict = {'AFR':'#f8766d','EUR':'#00b0f6','EAS':'#00bf7d','AMR':'#a3a500','SAS':'#e76bf3'}
		eur_pop_color_dict = {'TSI':'#e76bf3','FIN':'#a3a500','CEU':'#f8766d','GBR':'#00bf7d','IBS':'#00b0f6'}
		#giant - panel A
		labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
		indiv_pop = dict(zip(labeldata['sample'],labeldata['color']))
		labeldata = labeldata.rename(columns = {'sample':'IID'})
		fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
		df = pd.read_csv('../cache/inset_data/giant.height.eur.flashpcs',sep = '\t')
		df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
		eigenvals = pd.read_csv('../cache/inset_data/giant.height.eur.eigenvals',sep = '\t',header = None)
		eigenvals = (eigenvals/np.sum(eigenvals))*100
		ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
		ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
		ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
		plt.tight_layout()
		plt.savefig('../figures/main_text/giant.height.eur.PC1.v.PC2.pdf')
		plt.clf()
		#ukb fvc - panel b
		labeldata = pd.read_csv('../cache/inset_data/1kg_poplabel_map.txt', sep = '\t')
		labeldata['color'] = labeldata['pop'].map(eur_pop_color_dict)
		labeldata = labeldata.rename(columns = {'sample':'IID'})
		indiv_pop = dict(zip(labeldata['IID'],labeldata['color']))
		fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
		df = pd.read_csv('../cache/inset_data/ukb.fvc.eur.flashpcs',sep = '\t')
		df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
		eigenvals = pd.read_csv('../cache/inset_data/ukb.fvc.eur.eigenvals',sep = '\t',header = None)
		eigenvals = (eigenvals/np.sum(eigenvals))*100
		ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
		ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
		ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
		plt.tight_layout()
		plt.savefig('../figures/main_text/ukb.fvc.eur.PC1.v.PC2.pdf')
		plt.clf()
		#bolt-lmm household income - panel C
		labeldata = pd.read_csv('../cache/inset_data/1kg_poplabel_map.txt', sep = '\t')
		labeldata['color'] = labeldata['super_pop'].map(superpop_color_dict)
		labeldata = labeldata.rename(columns = {'sample':'IID'})
		indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
		fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
		df = pd.read_csv('../cache/inset_data/ukb.household_income.all.flashpcs',sep = '\t')
		df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
		eigenvals = pd.read_csv('../cache/inset_data/ukb.household_income.all.eigenvals',sep = '\t',header = None)
		eigenvals = (eigenvals/np.sum(eigenvals))*100
		ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
		ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
		ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
		plt.tight_layout()
		plt.savefig('../figures/main_text/ukb.household_income.all.PC1.v.PC2.pdf')
		plt.clf()
		#ukb waist circ - panel D
		labeldata = pd.read_csv('../cache/inset_data/1kg_poplabel_map.txt', sep = '\t')
		labeldata['color'] = labeldata['super_pop'].map(superpop_color_dict)
		labeldata = labeldata.rename(columns = {'sample':'IID'})
		indiv_superpop = dict(zip(labeldata['IID'],labeldata['color']))
		fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
		df = pd.read_csv('../cache/inset_data/ukb.waist_circ.all.flashpcs',sep = '\t')
		df = df.merge(labeldata[['IID','color']], on = 'IID', how = 'inner')
		eigenvals = pd.read_csv('../cache/inset_data/ukb.waist_circ.all.eigenvals',sep = '\t',header = None)
		eigenvals = (eigenvals/np.sum(eigenvals))*100
		ax.scatter(df['PC1'],df['PC2'], s=7, c = df['color'])
		ax.set_xlabel('PC1 (' + str(round(eigenvals.loc[0][0],2)) + '%)' )
		ax.set_ylabel('PC2 (' + str(round(eigenvals.loc[1][0],2)) + '%)')
		plt.tight_layout()
		plt.savefig('../figures/main_text/ukb.waist_circ.all.PC1.v.PC2.pdf')
		plt.clf()



	def plot_alphapillars(self):
		label = '1kg.all'
		alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
		alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
		alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
		alpha_dfplot = alpha_dfplot[['1e-05','1e-08','ycoordinate']]
		nsnp = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.nsnp.mat.txt',sep = '\t').set_index('Unnamed: 0')
		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (8,5), sharex = True)
		for i,j in enumerate(alpha_dfplot.columns[:2]):			
			ax[i].vlines(1, ymin = 0, ymax =17.5, color = 'black')
			ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
			ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
			ax[i].title.set_text('p < ' + str(j))
			ax[i].set_axisbelow(True)
			ax[i].set_xlabel('Isotropic inflation factor')
			ax[i].grid(True, linewidth = 1)
			if i == 0:
				ax[i].set_yticks(alpha_dfplot['ycoordinate'])
				ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
			else:
				ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
				ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
			ax[i].set_xticks([x for x in range(1,13)])
			ax[i].set_xticklabels([str(x) for x in ['1','','3','','5','','7','','9','','11','']])
		sns.despine()
		ax[0].set_ylim([0.5,17.5])
		ax[1].set_ylim([0.5,17.5])
		plt.tight_layout()
		plt.savefig('../figures/main_text/fig.wc.' + label + '.alpha.pillar.pdf')
		plt.clf()