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
import math
warnings.filterwarnings('ignore')

class stratification_comparison(object):

	def __init__(self, label_dict):

		self.label_dict = label_dict
		self.label_dict['mean'] = 'Mean value'
		self.label_dict['median'] = 'Median'
		self.label_dict[0] = 'Mean rank'

	def run(self):
		self.plot_comparison()
		
	def plot_comparison(self):
		for pc_count in ['6']:
			analyses = ['1kg.all']
			for analysis in analyses:
				trait_types ={'alcohol_intake_freq':'behavioral','household_income':'behavioral','neuroticism_score':'behavioral','overall_health':'behavioral','pack_years_smoking':'behavioral','years_schooling':'behavioral', \
				'birth_weight':'anthropometric','bmi':'anthropometric','height':'anthropometric','hip_circ':'anthropometric','waist_circ':'anthropometric', \
				'dbp':'other','fvc':'other','pulse_rate':'other','hand_grip_strength':'other','hair_color':'other','skin_color':'other'}

				wc_nopcs = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
				wc_ukb_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt', sep = '\t').set_index('Unnamed: 0')
				wc_1kg_pcs_only = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.' + analysis + '.pcs.only.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				wc_ukb_and_1kg_pcs = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.ukb.and.' + analysis + '.pcs.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				bolt_nopcs = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				bolt_wpcs = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')

				wc_nopcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.sad.pc.count.' + pc_count + '.mat.txt', sep = '\t').set_index('Unnamed: 0')
				wc_ukb_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.sad.pc.count.' + pc_count + '.mat.txt', sep = '\t').set_index('Unnamed: 0')
				wc_1kg_pcs_only_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.' + analysis + '.pcs.only.sps23.v2.sad.pc.count.' + pc_count + '.mat.txt',sep = '\t').set_index('Unnamed: 0')
				wc_ukb_and_1kg_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.ukb.and.' + analysis + '.pcs.sps23.v2.sad.pc.count.' + pc_count + '.mat.txt',sep = '\t').set_index('Unnamed: 0')
				bolt_nopcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.sad.pc.count.' + pc_count + '.mat.txt',sep = '\t').set_index('Unnamed: 0')
				bolt_wpcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.sad.pc.count.' + pc_count + '.mat.txt',sep = '\t').set_index('Unnamed: 0')

				def rand_jitter(arr,nanalysis,counter):	
					return arr - (0.9/nanalysis)*counter + 0.5

				methodlabels = ['No Adjustment','GWAS sample PCs','Prediction sample PCs','Both sets of PCs','Linear Mixed Model (LMM)','LMM + GWAS sample PCs']
				melt_methodlabels = {'wc.nopcs':6,'wc.pcs':5,'wc.1kg.pcs':4,'wc.ukb.and.1kg.pcs':3,'bolt.nopcs':2,'bolt.wpcs':1}
				for thresh in ['1e-05']:
					# fig, ax = plt.subplots(1,2, figsize = (10,5))

					# value = wc_nopcs_sads[[thresh]]
					# value.columns = ['wc.nopcs']
					# value['wc.pcs']=wc_ukb_pcs_sads.loc[value.index.tolist(),thresh]
					# value['wc.1kg.pcs']=wc_1kg_pcs_only_sads.loc[value.index.tolist(),thresh]
					# value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs_sads.loc[value.index.tolist(),thresh]
					# value['bolt.nopcs']=bolt_nopcs_sads.loc[value.index.tolist(),thresh]
					# value['bolt.wpcs']=bolt_wpcs_sads.loc[value.index.tolist(),thresh]
					# clip_value =value.clip(upper = 2.1,lower = 0.1)
					# melted = clip_value.melt()
					# melted['y'] = melted['variable'].map(melt_methodlabels) + np.random.uniform(-0.25,0.25,melted.shape[0])
					# ax[0].scatter(x=melted['value'],y=melted['y'], color = 'grey', alpha= 0.3, zorder = 0, s = 20)
					# means=value.mean(axis = 0)

					# ax[0].scatter(means, [x for x in range(len(value.columns),0,-1)], color = '#ca3a27', marker = '.', s = 200, zorder=3, label = '')
					# ax[0].set_yticks([x for x in range(len(value.columns),0,-1)])
					# ax[0].set_yticklabels(methodlabels,fontsize = 12)
					# ax[0].yaxis.set_tick_params("center")
					# ax[0].title.set_text('Number of significant PC-wise\nSAD components among top ' + pc_count + ' PCs')
					# ax[0].tick_params(axis='x', labelsize=12)

					# ax[0].set_xscale('log')
					# ax[0].set_xticks([0.1,0.25,0.5,0.75,1,2])
					# ax[0].set_xticklabels(['0','0.25','0.5','0.75','1','>2'])

					# value = wc_nopcs[[thresh]]
					# value.columns = ['wc.nopcs']
					# value['wc.pcs']=wc_ukb_pcs.loc[value.index.tolist(),thresh]
					# value['wc.1kg.pcs']=wc_1kg_pcs_only.loc[value.index.tolist(),thresh]
					# value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs.loc[value.index.tolist(),thresh]
					# value['bolt.nopcs']=bolt_nopcs.loc[value.index.tolist(),thresh]
					# value['bolt.wpcs']=bolt_wpcs.loc[value.index.tolist(),thresh]
					# clip_value =value.clip(upper = 5)
					# melted = clip_value.melt()
					# melted['y'] = melted['variable'].map(melt_methodlabels) + np.random.uniform(-0.25,0.25,melted.shape[0])

					# ax[1].scatter(x=melted['value'],y=melted['y'], color = 'grey', alpha= 0.3, zorder = 0, s = 20)
					# ax[1].scatter(value.mean(axis = 0).tolist(), [x for x in range(len(value.columns),0,-1)], color = '#CA6627', marker = '.', s = 200, zorder=3)					
					# ax[1].title.set_text('Isotropic inflation factor')
					# ax[1].set_yticklabels(['' for i in methodlabels],rotation=30,fontsize = 12)
					# ax[1].tick_params(axis='x', labelsize=12)

					# ax[1].set_xticks([1,2,3,4,5])
					# ax[1].set_xticklabels(['1','2','3','4','>5'])
					# sns.despine(top=True, right=True, left=False, bottom=False,)
					# plt.gcf().text(0.16,0.87, 'GWAS method to\nadjust for confounding', fontsize=14, ha = 'center')
					
					# plt.tight_layout()
					# plt.subplots_adjust(left=0.3)
					# plt.savefig('../figures/stratification_control_comp/strat.control.comparison.values.v2.pc.' + pc_count + '.' + analysis + '.pval.' + thresh + '.pdf')
					# plt.clf()


					#################make the same plot but using the ranks instead
					fig, ax = plt.subplots(1,2, figsize = (10,5), sharex = True)

					value = wc_nopcs_sads[[thresh]]
					value.columns = ['wc.nopcs']
					value['wc.pcs']=wc_ukb_pcs_sads.loc[value.index.tolist(),thresh]
					value['wc.1kg.pcs']=wc_1kg_pcs_only_sads.loc[value.index.tolist(),thresh]
					value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs_sads.loc[value.index.tolist(),thresh]
					value['bolt.nopcs']=bolt_nopcs_sads.loc[value.index.tolist(),thresh]
					value['bolt.wpcs']=bolt_wpcs_sads.loc[value.index.tolist(),thresh]
					value = value[['wc.nopcs','wc.pcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					
					value_rank = value.rank(axis=1, method = 'min')
					value_rank = (6-value_rank)+1
					value_rank = value_rank[['wc.nopcs','wc.pcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					rank_mean = pd.DataFrame(value_rank.mean(axis = 0)).T

					melted = value_rank.melt()
					melted['y'] = melted['variable'].map(melt_methodlabels) + np.random.uniform(-0.25,0.25,melted.shape[0])
					ax[0].scatter(x=melted['value'],y=melted['y'], color = 'grey', alpha= 0.3, zorder = 0, s = 20)
					print(6-rank_mean)
					ax[0].scatter(rank_mean, [x for x in range(len(value.columns),0,-1)], color = '#ca3a27', marker = '.', s = 200, zorder=3, label = '')
					ax[0].set_yticks([x for x in range(len(value.columns),0,-1)])
					ax[0].set_yticklabels(methodlabels,fontsize = 12)
					ax[0].yaxis.set_tick_params("center")
					ax[0].set_title('Significant PC-wise SAD components\namong top ' + pc_count + ' PCs', fontsize = 10)
					ax[0].tick_params(axis='x', labelsize=12)
					ax[0].set_xticks([1,2,3,4,5,6])
					ax[0].set_xticklabels(['6','5','4','3','2','1'])

					value = wc_nopcs[[thresh]]
					value.columns = ['wc.nopcs']
					value['wc.pcs']=wc_ukb_pcs.loc[value.index.tolist(),thresh]
					value['wc.1kg.pcs']=wc_1kg_pcs_only.loc[value.index.tolist(),thresh]
					value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs.loc[value.index.tolist(),thresh]
					value['bolt.nopcs']=bolt_nopcs.loc[value.index.tolist(),thresh]
					value['bolt.wpcs']=bolt_wpcs.loc[value.index.tolist(),thresh]
					value = value[['wc.nopcs','wc.pcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					value_rank = (6-value.rank(axis=1, method = 'min'))+1
					value_rank = value_rank[['wc.nopcs','wc.pcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					rank_mean = pd.DataFrame(value_rank.mean(axis = 0)).T
					print(6-rank_mean)
					melted = value_rank.melt()
					means=value.mean(axis = 0)
					###
					# print(means)
					melted['y'] = melted['variable'].map(melt_methodlabels) + np.random.uniform(-0.25,0.25,melted.shape[0])
					ax[1].scatter(x=melted['value'],y=melted['y'], color = 'grey', alpha= 0.3, zorder = 0, s = 20)
					ax[1].scatter(rank_mean, [x for x in range(len(value.columns),0,-1)], color = '#CA6627', marker = '.', s = 200, zorder=3)					
					ax[1].set_title('Isotropic inflation factor', fontsize = 10)
					ax[1].set_yticklabels(['' for i in methodlabels],rotation=30,fontsize = 12)
					ax[1].tick_params(axis='x', labelsize=12)
					ax[1].set_xticks([1,2,3,4,5,6])
					ax[1].set_xticklabels(['6','5','4','3','2','1'])

					sns.despine(top=True, right=True, left=False, bottom=False,)
					plt.gcf().text(0.16,0.87, 'GWAS method to\nadjust for confounding', fontsize=14, ha = 'center')
					fig.text(0.65, 0.01, 'Rank among methods', ha='center', fontsize =12)
					plt.tight_layout()
					plt.subplots_adjust(left=0.3)
					plt.savefig('../figures/stratification_control_comp/strat.control.comparison.rank.v2.pc.' + pc_count + '.' + analysis + '.pval.' + thresh + '.pdf')
					plt.clf()

					#heatmaps for supplement
					fig,ax = plt.subplots(2,3,figsize = (10,10), gridspec_kw = {'height_ratios' :[2/19.,17/19.],'width_ratios' :[8,0.5,8]})
					value = wc_ukb_pcs.sort_values(by=thresh)
					value = value[[thresh]]
					value.columns = ['wc.ukb']
					value['wc.nopcs']=wc_nopcs.loc[value.index.tolist(),thresh]
					value['wc.1kg.pcs']=wc_1kg_pcs_only.loc[value.index.tolist(),thresh]
					value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs.loc[value.index.tolist(),thresh]
					value['bolt.nopcs']=bolt_nopcs.loc[value.index.tolist(),thresh]
					value['bolt.wpcs']=bolt_wpcs.loc[value.index.tolist(),thresh]
					value['trait_type'] = value.index.to_series().map(trait_types)				
					value = value.drop(columns=['trait_type'])
					value_median = pd.DataFrame(value.median(axis=0),columns = ['median'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					value_mean = pd.DataFrame(value.mean(axis=0),columns = ['mean'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					value_groups = value_median.append(value_mean)
					value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					value = value.sort_values(by=['wc.ukb'])
					value = value.reindex(index=value.index[::-1])
					ordering =value.index.tolist()
					value = value.loc[ordering]
					value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					value_rank = value.rank(axis=1, method = 'min')-1
					value_rank = value_rank[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					rank_mean = pd.DataFrame(value_rank.mean(axis = 0)+1).T
					rank_mean = rank_mean.append(value_mean)
					rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
					rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]

					sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					ax[0,2].xaxis.set_label_position('top')
					ax[0,2].xaxis.tick_top()
					ax[0,2].set_xlabel('Isotropic inflation',fontsize = 14)
					ax[0,2].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
					ax[0,2].set_yticklabels(['' for trait in rank_mean.index.tolist()])
					ax[0,2].yaxis.set_tick_params(labelleft=False)
					ax[0,2].set_yticks([])
					ax[0,2].set_ylabel('')
					value_plot = value.rank(axis=1, method = 'min')-1
					sns.heatmap(value_plot,annot=False,ax=ax[1,2],cmap='Oranges',vmin = value.values.min(),vmax = value.values.max(),cbar = False)
					sns.heatmap(value_plot,annot=value, ax = ax[1,2],cmap='Oranges',fmt=".2f", cbar=False)
					ax[1,2].yaxis.set_tick_params(labelleft=False)
					ax[1,2].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
					ax[1,2].set_yticks([])
					ax[1,2].set_ylabel('')
					ax[1,2].set_xticks([])
					alphas = value
					ax[0,1].axis('off')

					for z,trait in enumerate(rank_mean_rank.index.tolist()[::-1]):
						ax[0,1].text(0.5,(z/2.)+0.25,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
					ax[1,1].axis('off')
					ax[1,1].set_xlim(ax[0,1].get_xlim())
					for z,trait in enumerate(value.index.tolist()[::-1]):
						ax[1,1].text(0.5,z/17.+0.03,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
					value = wc_ukb_pcs_sads.sort_values(by=thresh)
					value = value[[thresh]]
					value.columns = ['wc.ukb']
					value['wc.nopcs']=wc_nopcs_sads.loc[value.index.tolist(),thresh]
					value['wc.1kg.pcs']=wc_1kg_pcs_only_sads.loc[value.index.tolist(),thresh]
					value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs_sads.loc[value.index.tolist(),thresh]
					value['bolt.nopcs']=bolt_nopcs_sads.loc[value.index.tolist(),thresh]
					value['bolt.wpcs']=bolt_wpcs_sads.loc[value.index.tolist(),thresh]
					value['trait_type'] = value.index.to_series().map(trait_types)
					value = value.drop(columns=['trait_type'])
					value_mean = pd.DataFrame(value.mean(axis=0),columns = ['mean'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					value_median = pd.DataFrame(value.median(axis=0),columns = ['all'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					value_groups = value_median.append(value_mean)

					value_groups = value_groups.reindex(index=value_groups.index[::-1])
					value_groups = value_groups[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					value = value.loc[ordering]
					value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					value_rank = value.rank(axis=1, method = 'min')
					value_rank-=1
					rank_mean = pd.DataFrame(value_rank.mean(axis = 0)+1).T
					rank_mean = rank_mean.append(value_mean)
					rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
					rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]
					sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					ax[0,0].xaxis.set_label_position('top')
					ax[0,0].xaxis.tick_top()
					ax[0,0].set_xlabel('Number of significant PC-wise\nSAD components in top ' + pc_count + ' PCs',fontsize = 14)
					ax[0,0].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
					ax[0,0].set_yticklabels(['' for trait in rank_mean.index.tolist()])
					ax[0,0].yaxis.set_tick_params(labelleft=False)
					ax[0,0].set_yticks([])
					ax[0,0].set_ylabel('')
					value_plot = value.rank(axis=1, method = 'min')-1
					value = value.astype(int)

					sns.heatmap(value_plot,annot=False,ax=ax[1,0],cmap='Oranges',vmin = value.values.min(),vmax = value.values.max(),cbar = False)
					sns.heatmap(value_plot,annot=value, ax = ax[1,0],cmap='Oranges',fmt="d", cbar=False)
					ax[1,0].yaxis.set_tick_params(labelleft=False)
					ax[1,0].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
					ax[1,0].set_yticks([])
					ax[1,0].set_ylabel('')
					ax[1,0].set_xticks([])
					plt.suptitle('Comparison of SAD variance among methods to adjust for population structure',fontsize=16, y=1)

					plt.subplots_adjust(hspace=0.14)
					sns.despine(bottom=True, top=False, right=True,ax = ax[0,0])
					sns.despine(bottom=True, top=False, right=True,ax = ax[0,2])
					sns.despine(top=True, right=True, bottom=True, ax = ax[1,0])
					sns.despine(top=True, right=True, bottom=True, ax = ax[1,2])
					fig.subplots_adjust(top=0.8)
					plt.tight_layout()
					plt.savefig('../figures/stratification_control_comp/strat.control.comparison.heatmap.pc.' + pc_count + '.' + analysis + '.pval.' + thresh + '.pdf')
					plt.clf()

					# fig,ax = plt.subplots(2,3,figsize = (10,10), gridspec_kw = {'height_ratios' :[2/19.,17/19.],'width_ratios' :[4.5,1,4.5]})
					# value = wc_ukb_pcs.sort_values(by=thresh)
					# value = value[[thresh]]
					# value.columns = ['wc.ukb']
					# value['wc.nopcs']=wc_nopcs.loc[value.index.tolist(),thresh]
					# value['wc.1kg.pcs']=wc_1kg_pcs_only.loc[value.index.tolist(),thresh]
					# value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs.loc[value.index.tolist(),thresh]
					# value['bolt.nopcs']=bolt_nopcs.loc[value.index.tolist(),thresh]
					# value['bolt.wpcs']=bolt_wpcs.loc[value.index.tolist(),thresh]
					# value['trait_type'] = value.index.to_series().map(trait_types)				
					# value = value.drop(columns=['trait_type'])
					# value_median = pd.DataFrame(value.median(axis=0),columns = ['median'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					# value_mean = pd.DataFrame(value.mean(axis=0),columns = ['mean'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					# value_groups = value_median.append(value_mean)
					# value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# value = value.sort_values(by=['wc.ukb'])
					# value = value.reindex(index=value.index[::-1])
					# ordering =value.index.tolist()
					# value = value.loc[ordering]
					# value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# value_rank = value.rank(axis=1, method = 'min')-1
					# value_rank = value_rank[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# rank_mean = pd.DataFrame(value_rank.mean(axis = 0)+1).T
					# rank_mean = rank_mean.append(value_mean)
					# rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
					# rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]

					# sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					# sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,2],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					# ax[0,2].xaxis.set_label_position('top')
					# ax[0,2].xaxis.tick_top()
					# ax[0,2].set_xlabel('Isotropic inflation',fontsize = 14)
					# ax[0,2].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
					# ax[0,2].set_yticklabels(['' for trait in rank_mean.index.tolist()])
					# ax[0,2].yaxis.set_tick_params(labelleft=False)
					# ax[0,2].set_yticks([])
					# ax[0,2].set_ylabel('')
					# value_plot = value.rank(axis=1, method = 'min')-1
					# sns.heatmap(value_plot,annot=False,ax=ax[1,2],cmap='Oranges',vmin = value.values.min(),vmax = value.values.max(),cbar = False)
					# sns.heatmap(value_plot,annot=value, ax = ax[1,2],cmap='Oranges',fmt=".2f", cbar=False)
					# ax[1,2].yaxis.set_tick_params(labelleft=False)
					# ax[1,2].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
					# ax[1,2].set_yticks([])
					# ax[1,2].set_ylabel('')
					# ax[1,2].set_xticks([])
					# alphas = value
					# ax[0,1].axis('off')

					# for z,trait in enumerate(rank_mean_rank.index.tolist()[::-1]):
					# 	ax[0,1].text(0.5,(z/2.)+0.25,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
					# ax[1,1].axis('off')
					# ax[1,1].set_xlim(ax[0,1].get_xlim())
					# for z,trait in enumerate(value.index.tolist()[::-1]):
					# 	ax[1,1].text(0.5,z/17.+0.03,self.label_dict[trait], verticalalignment = 'center', horizontalalignment='center')
					# value = wc_ukb_pcs_sads.sort_values(by=thresh)
					# value = value[[thresh]]
					# value.columns = ['wc.ukb']
					# value['wc.nopcs']=wc_nopcs_sads.loc[value.index.tolist(),thresh]
					# value['wc.1kg.pcs']=wc_1kg_pcs_only_sads.loc[value.index.tolist(),thresh]
					# value['wc.ukb.and.1kg.pcs']=wc_ukb_and_1kg_pcs_sads.loc[value.index.tolist(),thresh]
					# value['bolt.nopcs']=bolt_nopcs_sads.loc[value.index.tolist(),thresh]
					# value['bolt.wpcs']=bolt_wpcs_sads.loc[value.index.tolist(),thresh]
					# value['trait_type'] = value.index.to_series().map(trait_types)
					# value = value.drop(columns=['trait_type'])
					# value_mean = pd.DataFrame(value.mean(axis=0),columns = ['mean'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					# value_median = pd.DataFrame(value.median(axis=0),columns = ['all'],index = ['wc.ukb','wc.nopcs','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']).T
					# value_groups = value_median.append(value_mean)

					# value_groups = value_groups.reindex(index=value_groups.index[::-1])
					# value_groups = value_groups[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# value = value.loc[ordering]
					# value = value[['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# value_rank = value.rank(axis=1, method = 'min')
					# value_rank-=1
					# rank_mean = pd.DataFrame(value_rank.mean(axis = 0)+1).T
					# rank_mean = rank_mean.append(value_mean)
					# rank_mean = rank_mean.loc[::-1][['wc.nopcs','wc.ukb','wc.1kg.pcs','wc.ukb.and.1kg.pcs','bolt.nopcs','bolt.wpcs']]
					# rank_mean_rank = rank_mean.rank(axis=1, method = 'min')-1
					# rank_mean_rank.loc['mean'] = rank_mean_rank.loc[0]
					# sns.heatmap(rank_mean_rank,annot=False,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					# sns.heatmap(rank_mean_rank,annot=rank_mean,ax=ax[0,0],cmap='Oranges',vmin = rank_mean_rank.values.min(),vmax = rank_mean_rank.values.max(),cbar = False)
					# ax[0,0].xaxis.set_label_position('top')
					# ax[0,0].xaxis.tick_top()
					# ax[0,0].set_xlabel('Number of significant PC-wise\nSAD components in top ' + pc_count + ' PCs',fontsize = 14)
					# ax[0,0].set_xticklabels(methodlabels,rotation=30,ha='left',fontsize = 9)
					# ax[0,0].set_yticklabels(['' for trait in rank_mean.index.tolist()])
					# ax[0,0].yaxis.set_tick_params(labelleft=False)
					# ax[0,0].set_yticks([])
					# ax[0,0].set_ylabel('')
					# value_plot = value.rank(axis=1, method = 'min')-1
					# value = value.astype(int)

					# sns.heatmap(value_plot,annot=False,ax=ax[1,0],cmap='Oranges',vmin = value.values.min(),vmax = value.values.max(),cbar = False)
					# sns.heatmap(value_plot,annot=value, ax = ax[1,0],cmap='Oranges',fmt="d", cbar=False)
					# ax[1,0].yaxis.set_tick_params(labelleft=False)
					# ax[1,0].xaxis.set_tick_params(labelbottom=False, labeltop = False, top = False, bottom=False)
					# ax[1,0].set_yticks([])
					# ax[1,0].set_ylabel('')
					# ax[1,0].set_xticks([])
					# plt.suptitle('Comparison of SAD variance among methods to adjust for population structure',fontsize=16, y=1)
					# plt.tight_layout()
					# plt.subplots_adjust(hspace=0.14)
					# sns.despine(bottom=True, top=False, right=True,ax = ax[0,0])
					# sns.despine(bottom=True, top=False, right=True,ax = ax[0,2])
					# sns.despine(top=True, right=True, bottom=True, ax = ax[1,0])
					# sns.despine(top=True, right=True, bottom=True, ax = ax[1,2])
					# fig.subplots_adjust(top=0.8)
					# plt.savefig('../figures/stratification_control_comp/strat.control.comparison.pc.' + pc_count + '.' + analysis + '.pval.' + thresh + '.pdf')
					# plt.clf()
