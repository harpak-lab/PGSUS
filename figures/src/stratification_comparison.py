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

class stratification_comparison(object):

	def __init__(self, analyses, label_dict):

		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		self.plot_comparison(self)
		
	def plot_comparison(self):
		analyses = ['1kg.all','1kg.eur']
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

			wc_nopcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt', sep = '\t').set_index('Unnamed: 0')
			wc_ukb_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.sps23.aperm.1K.to.1M.v2.sad.pc.count.mat.txt', sep = '\t').set_index('Unnamed: 0')
			wc_1kg_pcs_only_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.' + analysis + '.pcs.only.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
			wc_ukb_and_1kg_pcs_sads = pd.read_csv('../cache/alpha_matrices/plink.wc.' + analysis + '.ukb.and.' + analysis + '.pcs.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
			bolt_nopcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')
			bolt_wpcs_sads = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + analysis + '.sps23.v2.sad.pc.count.mat.txt',sep = '\t').set_index('Unnamed: 0')

			def rand_jitter(arr,nanalysis,counter):	
				return arr - (0.9/nanalysis)*counter + 0.5

			methodlabels = ['None','GWAS sample PCs','Target sample PCs','Both sets of PCs','Linear Mixed Model\n(LMM)','LMM with GWAS\nsample PCs']
			for thresh in ['1.0','0.001','1e-05','1e-08']:
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
