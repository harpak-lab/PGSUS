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

class decomps_nondirect(object):

	def __init__(self, analyses, label_dict, sps_traits):

		self.analyses = analyses
		self.label_dict = label_dict
		self.sps_traits = sps_traits
		self.palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41', 'nondirect':'#8a461b'}
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
		if 'giant2022' in self.analyses:
			self.plot_giant2022()
		if 'okbay2022' in self.analyses:
			self.plot_okbay2022()
		if 'akbari2024' in self.analyses:
			self.plot_akbari2024()

	def plot_akbari2024(self):
		row_counter = 0
		fig, ax = plt.subplots(nrows = 4, ncols = 4, figsize = (4*(self.pcs*0.8+1),20))
		for trait in ['household_income', 'overall_health', 'skin_color', 'years_schooling']:
			column_counter = 0
			for pval in self.pval_array:
				col_index = column_counter % 4
				df = pd.read_csv('../cache/component_inputs/nondirect/akbari2024/akbari2024.' + trait + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t').set_index('Unnamed: 0')
				pcs=self.pcs
				palette = self.palette
				col_nums = np.array([int(i+1) for i in range(pcs)])

				ax[row_counter,col_index].fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
				ax[row_counter,col_index].fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
				ax[row_counter,col_index].fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
				ax[row_counter,col_index].fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

				ax[row_counter,col_index].scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
				ax[row_counter,col_index].scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
				ax[row_counter,col_index].scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
				ax[row_counter,col_index].scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
				
				df = df[df.columns[:6]]
				df.columns = [i for i in range(6)]
				direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
				sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
				covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
				nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

				if len(direct_sig_cols) != 0:
					ax[row_counter,col_index].scatter(col_nums[direct_sig_cols]-0.25, df.loc['direct_vc_estimate'][direct_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['direct'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
				if len(sad_sig_cols) != 0:
					ax[row_counter,col_index].scatter(col_nums[sad_sig_cols]-0.083, df.loc['sad_vc_estimate'][sad_sig_cols].tolist(), marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['sad'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
				if len(covar_sig_cols) != 0:
					ax[row_counter,col_index].scatter(col_nums[covar_sig_cols]+0.083, df.loc['covar_vc_estimate'][covar_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['covar'], linewidth = 2, linestyle = (0, (1, 1)), label = '')
				if len(nondirect_sig_cols) != 0:
					ax[row_counter,col_index].scatter(col_nums[nondirect_sig_cols]+0.25, df.loc['nondirect_vc_estimate'][nondirect_sig_cols], marker = 'o', s = 200, facecolor = 'none', edgecolor = palette['nondirect'], linewidth = 2, linestyle = (0, (1, 1)), label = '')

				upper_y = ax[row_counter,col_index].get_ylim()[1]
				lower_y = ax[row_counter,col_index].get_ylim()[0]
				asterisk_y = np.array([upper_y for i in range(df.shape[1])])
				ax[row_counter,col_index].hlines(0, 0, np.max(col_nums)+1, 'grey','-', zorder = 0)
				for i in col_nums:
					ax[row_counter,col_index].vlines(i+0.5, lower_y, upper_y, 'grey', 'dotted', zorder = 0)
				ax[row_counter,col_index].set_xticks([int(i+1) for i in range(pcs)])
				ax[row_counter,col_index].set_xticklabels(['PC' + str(i+1) for i in range(pcs)])
				ax[row_counter,col_index].set_xlim(0.5, pcs+0.5)
				ax[row_counter,col_index].set_ylabel('Component / Total PGS Variance')
				ax[row_counter,col_index].title.set_text(r'$p$-value $<$ ' + str(pval))
				column_counter += 1
			row_counter +=1
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/component_plots_sps_v2/nondirect/akbari2024/combined.component.decomp.plots.' + str(self.pcs) + '.pcs.pdf')
		plt.clf()

	def plot_okbay2022(self):
		for pval in self.pval_array:
			for label in ['1kg.all', '1kg.eur']:
				df = pd.read_csv('../cache/component_inputs/nondirect/okbay2022/okbay.' + label + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t').set_index('Unnamed: 0')
				fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
				pcs=self.pcs
				palette = self.palette
				col_nums = np.array([int(i+1) for i in range(pcs)])

				ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
				ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
				ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
				ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')
				print(df)
				ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
				ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
				ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
				ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
				
				df = df[df.columns[:6]]
				df.columns = [i for i in range(6)]
				direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
				sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
				covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
				nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
				ax.set_ylabel('Component / Total PGS Variance')

				plt.legend(loc='lower right')
				plt.savefig('../figures/component_plots_sps_v2/nondirect/okbay2022/component.plot.' + label +'.pval.' + str(pval) + '.' + str(self.pcs) + '.pcs.pdf')
				plt.clf()

	def plot_giant2022(self):
		for training in ['all','eur']:
			for pval in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:

					df = pd.read_csv('../cache/component_inputs/nondirect/giant2022/giant.' + training + '.' + label + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/giant2022/component.plot.' + training + '.' + label +'.pval.' + str(pval) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_giant(self):
		for trait in ['giant_height','giant_height_rescaled']:
			for pval in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:

					df = pd.read_csv('../cache/component_inputs/nondirect/giant/' + trait + '.' + label + '.block.permutation.stats.pval.' + str(pval) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/giant/component.plot.' + label + '.' + trait +'.pval.' + str(pval) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_ukb_and_1kg(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for cohort,label in zip(['1kg.all', '1kg.eur'],['ukb.and.1kg.all.pcs','ukb.and.1kg.eur.pcs']):
					df = pd.read_csv('../cache/component_inputs/nondirect/ukb.and.1kg.pcs/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')

					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/ukb.and.1kg.pcs/component.plot.' + cohort + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_1kg_pcs_only(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for cohort,label in zip(['1kg.all', '1kg.eur'],['1kg.all.pcs.only','1kg.eur.pcs.only']):
					df = pd.read_csv('../cache/component_inputs/nondirect/1kg.pcs.only/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/1kg.pcs.only/component.plot.' + cohort + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_bolt(self):
		for grm in ['wpcs','nopcs']:
			for trait in self.sps_traits:
				for thresh in self.pval_array:
					for label in ['1kg.all', '1kg.eur']:
						df = pd.read_csv('../cache/component_inputs/nondirect/bolt/bolt.' + grm + '.' + label + '.' + trait + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
						fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
						pcs=self.pcs
						palette = self.palette
						col_nums = np.array([int(i+1) for i in range(pcs)])

						ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
						ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

						ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
						ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
						ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
						ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
						
						df = df[df.columns[:6]]
						df.columns = [i for i in range(6)]
						direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
						sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
						covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
						nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
						ax.set_ylabel('Component / Total PGS Variance')

						plt.legend(loc='lower right')
						plt.savefig('../figures/component_plots_sps_v2/nondirect/bolt/component.plot.bolt.' + grm + '.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
						plt.clf()


	def plot_nopcs(self):
		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.all', '1kg.eur']:

					df = pd.read_csv('../cache/component_inputs/nondirect/nopcs/plink.wc.nopcs.' + label + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])

					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] > df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]

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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/nopcs/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

	def plot_wc(self):
		counter = 0
		total_negative_covar = 0
		total_negative_nondirect = 0
		total_positive_sad = 0
		total_pos_sad_neg_nondirect = 0
		total_pos_sad_pos_nondirect = 0
		total_neg_covar_neg_nondirect = 0
		total_neg_covar_pos_nondirect = 0

		#significant counts
		sig_total_negative_covar = 0
		sig_total_negative_nondirect = 0
		sig_total_positive_sad = 0
		sig_total_pos_sad_neg_nondirect = 0
		sig_total_pos_sad_pos_nondirect = 0
		sig_total_neg_covar_neg_nondirect = 0
		sig_total_neg_covar_pos_nondirect = 0

		for trait in self.sps_traits:
			for thresh in self.pval_array:
				for label in ['1kg.eur','1kg.all']:
					df = pd.read_csv('../cache/component_inputs/nondirect/wc/plink.wc.' + label + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (self.pcs*0.8+1,5))
					pcs=self.pcs
					palette = self.palette
					col_nums = np.array([int(i+1) for i in range(pcs)])
					# print(np.where(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]))
					ax.fill_between(col_nums-0.25, df.loc['upper95_perm_direct'].tolist()[:pcs], df.loc['lower0_perm_direct'].tolist()[:pcs], where=(df.loc['upper95_perm_direct'].tolist()[:pcs] > df.loc['lower0_perm_direct'].tolist()[:pcs]), facecolor=palette['direct'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums-0.083, df.loc['upper95_perm_sad'].tolist()[:pcs], df.loc['lower0_perm_sad'].tolist()[:pcs], where=(df.loc['upper95_perm_sad'].tolist()[:pcs] > df.loc['lower0_perm_sad'].tolist()[:pcs]), facecolor=palette['sad'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.083, df.loc['upper975_perm_covar'].tolist()[:pcs], df.loc['lower025_perm_covar'].tolist()[:pcs], where=(df.loc['upper975_perm_covar'].tolist()[:pcs] > df.loc['lower025_perm_covar'].tolist()[:pcs]), facecolor=palette['covar'], alpha=0.2, edgecolor = 'none')
					ax.fill_between(col_nums+0.25, df.loc['upper975_perm_nondirect'].tolist()[:pcs], df.loc['lower025_perm_nondirect'].tolist()[:pcs], where=(df.loc['upper975_perm_nondirect'].tolist()[:pcs] < df.loc['lower025_perm_nondirect'].tolist()[:pcs]), facecolor=palette['nondirect'], alpha=0.2, edgecolor = 'none')

					ax.scatter(col_nums-0.25, df.loc['direct_vc_estimate'].tolist()[:pcs], color = palette['direct'], label='direct variance')
					ax.scatter(col_nums-0.083, df.loc['sad_vc_estimate'].tolist()[:pcs], color = palette['sad'], label='SAD variance')
					ax.scatter(col_nums+0.083, df.loc['covar_vc_estimate'].tolist()[:pcs], color = palette['covar'], label='direct-SAD covariance')
					ax.scatter(col_nums+0.25, df.loc['nondirect_vc_estimate'].tolist()[:pcs], color = palette['nondirect'], label='nondirect variance')
					
					df = df[df.columns[:6]]
					df.columns = [i for i in range(6)]
					direct_sig_cols = df.columns[df.loc['direct_vc_pvals'][:6] < 0.05]
					sad_sig_cols = df.columns[df.loc['sad_vc_pvals'][:6] < 0.05]
					covar_sig_cols = df.columns[df.loc['covar_vc_pvals'][:6] < 0.025]
					nondirect_sig_cols = df.columns[df.loc['nondirect_vc_pvals'][:6] < 0.025]
					
					sad_temp = df.columns[(df.loc['sad_vc_pvals'][:6] < 0.05)& (df.loc['sad_vc_estimate'] > 0)]
					covar_temp = df.columns[(df.loc['covar_vc_pvals'][:6] < 0.025)& (df.loc['covar_vc_estimate'] < 0)]
					nondirect_temp = df.columns[(df.loc['nondirect_vc_pvals'][:6] < 0.025)& (df.loc['nondirect_vc_estimate'] < 0)]
					sig_total_negative_covar += covar_temp.shape[0]
					sig_total_negative_nondirect += nondirect_temp.shape[0]

					nondirect_pos_temp = df.columns[(df.loc['nondirect_vc_pvals'][:6] < 0.025)& (df.loc['nondirect_vc_estimate'] > 0)]
					nondirect_neg_temp = df.columns[(df.loc['nondirect_vc_pvals'][:6] < 0.025)& (df.loc['nondirect_vc_estimate'] < 0)]
					

					sig_total_positive_sad += sad_temp.shape[0]
					
					sad_nondirect_pos = np.in1d(sad_temp,nondirect_pos_temp).sum()
					sad_nondirect_neg = np.in1d(sad_temp,nondirect_neg_temp).sum()
					sig_total_pos_sad_pos_nondirect += np.in1d(sad_temp,nondirect_pos_temp).sum()
					sig_total_pos_sad_neg_nondirect += np.in1d(sad_temp,nondirect_neg_temp).sum()
					sig_total_neg_covar_neg_nondirect += np.in1d(covar_temp,nondirect_neg_temp).sum()
					sig_total_neg_covar_pos_nondirect += np.in1d(covar_temp,nondirect_pos_temp).sum()

					#add to total number depending
					total_positive_sad += (df.loc['sad_vc_estimate'] > 0).sum()
					total_negative_covar += (df.loc['covar_vc_estimate'] < 0).sum()
					total_negative_nondirect += (df.loc['nondirect_vc_estimate'] < 0).sum()
					total_pos_sad_neg_nondirect += np.in1d(df.columns[(df.loc['nondirect_vc_estimate'] <= 0)],df.columns[df.loc['sad_vc_estimate'] > 0]).sum()
					total_pos_sad_pos_nondirect += np.in1d(df.columns[(df.loc['nondirect_vc_estimate'] > 0)],df.columns[df.loc['sad_vc_estimate'] > 0]).sum()
					total_neg_covar_neg_nondirect += np.in1d(df.columns[(df.loc['covar_vc_estimate'] <= 0)],df.columns[df.loc['nondirect_vc_estimate'] < 0]).sum()
					total_neg_covar_pos_nondirect += np.in1d(df.columns[(df.loc['covar_vc_estimate'] <= 0)],df.columns[df.loc['nondirect_vc_estimate'] > 0]).sum()
					#add to significant totals


					counter +=1


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
					ax.set_ylabel('Component / Total PGS Variance')

					plt.legend(loc='lower right')
					plt.savefig('../figures/component_plots_sps_v2/nondirect/wc/component.plot.' + label + '.' + trait + '.pval.' + str(thresh) + '.' + str(self.pcs) + '.pcs.pdf')
					plt.clf()

		# print(total_positive_sad, total_negative_covar, total_negative_nondirect, total_pos_sad_neg_nondirect, total_pos_sad_pos_nondirect, total_neg_covar_neg_nondirect, total_neg_covar_pos_nondirect)
		# print(counter)
		# print(sig_total_positive_sad, sig_total_pos_sad_neg_nondirect, sig_total_pos_sad_pos_nondirect, sig_total_negative_covar, sig_total_negative_nondirect, sig_total_neg_covar_neg_nondirect,sig_total_neg_covar_pos_nondirect)











