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

class decomps(object):

	def __init__(self, analyses, label_dict, sps_traits):

		self.analyses = analyses
		self.label_dict = label_dict
		self.sps_traits = sps_traits
		self.palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41'}
		self.pcs = 6
		self.pval_array = [1.0,0.001,0.00001,0.00000001]
	
	def run(self):
		print(self.analyses)
		if 'wc' in self.analyses:
			self.plot_wc()
		if 'nopcs' in self.analyses:
			self.plot_nopcs()
		if 'bolt' in self.analyses:
			self.plot_bolt()
		if '1kg.pcs.only' in self.analyses:
			self.plot_1kg_pcs_only()
		if 'ukb.and.1kg.pcs' in self.analyses:
			self.plot_ukb_and_1kg()
		if 'giant' in self.analyses:
			self.plot_giant()
		if 'maf01' in self.analyses:
			self.plot_maf01()
		if 'genotypes' in self.analyses:
			self.plot_genotypes()
	
	def plot_giant(self):
		for trait in ['giant_height','giant_height_rescaled']:
			for pval in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:
					print(trait,pval,label)
					df = pd.read_csv('../cache/component_inputs/giant/' + trait + '/giant.' + label + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])
					df = df.set_index('Unnamed: 0')

					ax.fill_between(col_nums-0.25, df.loc['upper975_perm_direct'].to_numpy()[:self.pcs], df.loc['lower025_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower025_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper975_perm_sad'].to_numpy()[:self.pcs], df.loc['lower025_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower025_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/giant/component.plot.' + label + '.' + trait +'.pval.' + str(pval) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_ukb_and_1kg(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for cohort,label in zip(['1kg.all', '1kg.eur'],['ukb.and.1kg.all.pcs','ukb.and.1kg.eur.pcs']):
					df = pd.read_csv('../cache/component_inputs/ukb.and.1kg.pcs/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					print(label,cohort,trait,thresh)
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/ukb.and.1kg.pcs/component.plot.' + cohort + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_1kg_pcs_only(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for cohort,label in zip(['1kg.all', '1kg.eur'],['1kg.all.pcs.only','1kg.eur.pcs.only']):
					df = pd.read_csv('../cache/component_inputs/1kg.pcs.only/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					print(label,cohort,trait,thresh)
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/1kg.pcs.only/component.plot.' + cohort + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_bolt(self):
		for grm in ['wpcs','nopcs']:
			for trait in self.sps_traits:
				for thresh in self.pval_array:
					for label in ['1kg.all', '1kg.eur']:
						df = pd.read_csv('../cache/component_inputs/bolt/bolt.' + grm + '.' + label + '.' + trait + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
						print('bolt.'+grm,trait,thresh,label)
						fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
						col_nums = np.array([int(i+1) for i in range(self.pcs)])

						ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
						
						ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
						ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
						ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
						
						upper_y = ax.get_ylim()[1]
						lower_y = ax.get_ylim()[0]
						asterisk_y = np.array([upper_y for i in range(df.shape[1])])

						ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
						for i in col_nums:
							ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
						
						ax.set_xticks([int(i+1) for i in range(self.pcs)])
						ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
						ax.set_xlim(0.5, self.pcs+0.5)
						ax.set_ylabel('Component / Total PGS Variance')
						sns.despine()

						plt.legend(loc='lower right')
						plt.savefig('../figures/component_plots_sps_v2/bolt/component.plot.bolt.' + grm + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
						plt.clf()

	def plot_halves(self):
		for half in ['1','2']:
			for trait in self.sps_traits:
				for thresh in self.pval_array:
					for label in ['1kg.all', '1kg.eur']:
						df = pd.read_csv('../cache/component_inputs/ascertain_validate/plink.half.' + half + '.' + label + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
						print('half' + half,trait,thresh,label)
						fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
						col_nums = np.array([int(i+1) for i in range(self.pcs)])
						ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
						ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
						ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
						ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
						upper_y = ax.get_ylim()[1]
						lower_y = ax.get_ylim()[0]
						asterisk_y = np.array([upper_y for i in range(df.shape[1])])
						ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
						for i in col_nums:
							ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
						ax.set_xticks([int(i+1) for i in range(self.pcs)])
						ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
						ax.set_xlim(0.5, self.pcs+0.5)
						ax.set_ylabel('Component / Total PGS Variance')
						sns.despine()
						plt.legend(loc='lower right')
						plt.savefig('../figures/component_plots_sps_v2/ascertain_validate/component.plot.half.' + half + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
						plt.clf()

	def plot_nosingletons(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:
					print('nosingletons',trait,thresh,label)
					df = pd.read_csv('../cache/component_inputs/nosingletons/plink.wc.' + label + '.' + trait + '.nosingletons.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nosingletons/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_nopcs(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:
					print('nopcs',trait,thresh,label)
					df = pd.read_csv('../cache/component_inputs/nopcs/plink.wc.nopcs.' + label + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])
					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])
					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()
					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nopcs/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_wc(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:
					print('wc',trait,thresh,label)
					df = pd.read_csv('../cache/component_inputs/wc/plink.wc.' + label + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/wc/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_maf01(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.eur','1kg.all']:
					print('maf01',trait,thresh,label)
					df = pd.read_csv('../cache/component_inputs/maf01/plink.wc.' + label + '.sps23.' + trait + '.maf01.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])
					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/maf01/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_genotypes(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.eur','1kg.all']:
					print('genotypes',trait,thresh,label)
					df = pd.read_csv('../cache/component_inputs/genotypes/plink.wc.' + label + '.sps23.' + trait + '.genotypes.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.67+1,5))
					col_nums = np.array([int(i+1) for i in range(self.pcs)])
					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].to_numpy()[:self.pcs], df.loc['lower0_perm_direct'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_direct'].to_numpy()[:self.pcs] > df.loc['lower0_perm_direct'].to_numpy()[:self.pcs]), facecolor=self.palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_covar'].to_numpy()[:self.pcs], df.loc['lower025_perm_covar'].to_numpy()[:self.pcs], where=(df.loc['upper975_perm_covar'].to_numpy()[:self.pcs] > df.loc['lower025_perm_covar'].to_numpy()[:self.pcs]), facecolor=self.palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums, df.loc['upper95_perm_sad'].to_numpy()[:self.pcs], df.loc['lower0_perm_sad'].to_numpy()[:self.pcs], where=(df.loc['upper95_perm_sad'].to_numpy()[:self.pcs] > df.loc['lower0_perm_sad'].to_numpy()[:self.pcs]), facecolor=self.palette['sad'], alpha=0.2, edgecolor = 'none')
					
					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['direct'], label='direct variance')
					ax.scatter(col_nums, df.loc['sad_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.25, df.loc['covar_vc_estimate'].to_numpy()[:self.pcs], color = self.palette['covar'], label='direct-SAD covariance')
					
					upper_y = ax.get_ylim()[1]
					lower_y = ax.get_ylim()[0]
					asterisk_y = np.array([upper_y for i in range(df.shape[1])])

					ax.hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
					for i in col_nums:
						ax.vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
					
					ax.set_xticks([int(i+1) for i in range(self.pcs)])
					ax.set_xticklabels(['PC' + str(i+1) for i in range(self.pcs)])
					ax.set_xlim(0.5, self.pcs+0.5)
					ax.set_ylabel('Component / Total PGS Variance')
					sns.despine()

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/genotypes/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()













