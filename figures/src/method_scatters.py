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

class method_scatters(object):

	def __init__(self, analyses, label_dict):
		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		if 'wc_v_nopcs' in analyses:
			self.plot_wc_v_nopcs()
		if 'wc_v_nosingletons' in analyses:
			self.plot_wc_v_nosingletons()
		if 'wc_v_bolt' in analyses:
			self.plot_wc_v_bolt()
		if 'bolt_v_bolt' in analyses:
			self.plot_bolt_v_bolt()
		if 'ascertainment_v_validation' in analyses:
			self.plot_ascertainment_v_validation()
		if 'wc_v_shuffle' in analyses:
			self.plot_wc_v_shuffle()
	
	def estimate_gaussian_outliers(dataset):
		mu = np.mean(dataset)
		sigma = np.std(dataset)
		limit = sigma * 1.5
		min_threshold = mu - limit
		max_threshold = mu + limit
		return mu, sigma, min_threshold, max_threshold

	def plot_wc_v_nopcs(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_nopcs = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_nopcs = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_nopcs = alpha_df_nopcs.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_nopcs = alpha_se_df_nopcs.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot[j].values)
				condition_alpha_dfplot = (alpha_dfplot[j] > max_threshold)
				outliers_alpha_dfplot = alpha_dfplot[j].loc[condition_alpha_dfplot]
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot_nopcs[j].values)
				condition_alpha_dfplot_nopcs = (alpha_dfplot_nopcs[j] > max_threshold)
				outliers_alpha_dfplot_nopcs = alpha_dfplot_nopcs[j].loc[condition_alpha_dfplot_nopcs]
				all_outliers = set(outliers_alpha_dfplot.index.tolist() + outliers_alpha_dfplot_nopcs.index.tolist())
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_nopcs[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_nopcs[j].loc[all_outliers], color = '#CA6627', s = 100, marker = '*')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_nopcs[j], xerr = alpha_se_dfplot[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_nopcs[j], yerr = alpha_se_dfplot_nopcs[j], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$ with PCs')
				ax[i].set_ylabel(r'$\hat{\alpha}$ without PCs')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if alpha_dfplot[j].loc[outlier] > alpha_dfplot_nopcs[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_nopcs[j].loc[outlier],label_dict[outlier],color = '#CA6627',horizontalalignment='right',verticalalignment='top')
					if alpha_dfplot[j].loc[outlier] < alpha_dfplot_nopcs[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_nopcs[j].loc[outlier],label_dict[outlier],color = '#CA6627')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/method_scatters/wc.v.nopcs.' + label + '.alpha.scatter.pdf')
			plt.clf()

	def plot_wc_v_shuffle(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_shuffle = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_shuffle = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_shuffle = alpha_df_shuffle.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_shuffle = alpha_se_df_shuffle.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot[j].values)
				condition_alpha_dfplot = (alpha_dfplot[j] > max_threshold)
				outliers_alpha_dfplot = alpha_dfplot[j].loc[condition_alpha_dfplot]
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot_shuffle[j].values)
				condition_alpha_dfplot_shuffle = (alpha_dfplot_shuffle[j] > max_threshold)
				outliers_alpha_dfplot_shuffle = alpha_dfplot_shuffle[j].loc[condition_alpha_dfplot_shuffle]
				all_outliers = set(outliers_alpha_dfplot.index.tolist() + outliers_alpha_dfplot_shuffle.index.tolist())
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_shuffle[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_shuffle[j].loc[all_outliers], color = '#CA6627', s = 100, marker = '*')				
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_shuffle[j], xerr = alpha_se_dfplot[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_shuffle[j], yerr = alpha_se_dfplot_shuffle[j], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$ with PCs')
				ax[i].set_ylabel(r'$\hat{\alpha}$ with shuffled $p$-values')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if alpha_dfplot[j].loc[outlier] > alpha_dfplot_shuffle[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_shuffle[j].loc[outlier],label_dict[outlier],color = '#CA6627',horizontalalignment='right',verticalalignment='top')
					if alpha_dfplot[j].loc[outlier] < alpha_dfplot_shuffle[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_shuffle[j].loc[outlier],label_dict[outlier],color = '#CA6627')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/method_scatters/wc.v.shuffle.' + label + '.alpha.scatter.pdf')
			plt.clf()

	def plot_wc_v_nosingletons(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_nosingletons = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.nosingletons.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_nosingletons = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.nosingletons.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_nosingletons = alpha_df_nosingletons.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_nosingletons = alpha_se_df_nosingletons.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot[j].values)
				condition_alpha_dfplot = (alpha_dfplot[j] > max_threshold)
				outliers_alpha_dfplot = alpha_dfplot[j].loc[condition_alpha_dfplot]
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot_nosingletons[j].values)
				condition_alpha_dfplot_nosingletons = (alpha_dfplot_nosingletons[j] > max_threshold)
				outliers_alpha_dfplot_nosingletons = alpha_dfplot_nosingletons[j].loc[condition_alpha_dfplot_nosingletons]
				all_outliers = set(outliers_alpha_dfplot.index.tolist() + outliers_alpha_dfplot_nosingletons.index.tolist())
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_nosingletons[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_nosingletons[j].loc[all_outliers], color = '#CA6627', s = 100, marker = '*')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_nosingletons[j], xerr = alpha_se_dfplot[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_nosingletons[j], yerr = alpha_se_dfplot_nosingletons[j], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$ with all clumps')
				ax[i].set_ylabel(r'$\hat{\alpha}$ without single SNP clumps')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if alpha_dfplot[j].loc[outlier] > alpha_dfplot_nosingletons[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_nosingletons[j].loc[outlier],label_dict[outlier],color = '#CA6627',horizontalalignment='right',verticalalignment='top')
					if alpha_dfplot[j].loc[outlier] < alpha_dfplot_nosingletons[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_nosingletons[j].loc[outlier],label_dict[outlier],color = '#CA6627')


			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/comparison_scatters/wc.v.nosingletons.' + label + '.alpha.scatter.pdf')
			plt.clf()

	def plot_wc_v_bolt(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_bolt = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_bolt = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_bolt = alpha_df_bolt.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_bolt = alpha_se_df_bolt.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				all_outliers = alpha_dfplot.loc[(alpha_dfplot_bolt[j] > alpha_dfplot[j])].index.tolist()
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_bolt[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_bolt[j].loc[all_outliers], color = 'black', s = 100, marker = '*')
				ax[i].errorbar(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_bolt[j].loc[~alpha_dfplot.index.isin(all_outliers)], xerr = alpha_se_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_bolt[j].loc[~alpha_dfplot.index.isin(all_outliers)], yerr = alpha_se_dfplot_bolt[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$')
				ax[i].set_ylabel(r'$\hat{\alpha}$ BOLT-LMM')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if outlier == 'years_schooling':
						if alpha_dfplot[j].loc[outlier] < alpha_dfplot_bolt[j].loc[outlier]:
							ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_bolt[j].loc[outlier]-0.2,label_dict[outlier],color = 'black',verticalalignment='top')
					else:
						if alpha_dfplot[j].loc[outlier] > alpha_dfplot_bolt[j].loc[outlier]:
							ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_bolt[j].loc[outlier],label_dict[outlier],color = 'black',horizontalalignment='right',verticalalignment='top')
						if alpha_dfplot[j].loc[outlier] < alpha_dfplot_bolt[j].loc[outlier]:
							ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_bolt[j].loc[outlier],label_dict[outlier],color = 'black',verticalalignment='top')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/comparison_scatters/wc.v.bolt.' + label + '.alpha.scatter.pdf')
			plt.clf()

	def plot_bolt_v_bolt(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/bolt.nopcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_bolt = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_bolt = pd.read_csv('../cache/alpha_matrices/bolt.wpcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_bolt = alpha_df_bolt.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_bolt = alpha_se_df_bolt.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot[j].values)
				condition_alpha_dfplot = (alpha_dfplot[j] > max_threshold)
				outliers_alpha_dfplot = alpha_dfplot[j].loc[condition_alpha_dfplot]
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot_bolt[j].values)
				condition_alpha_dfplot_bolt = (alpha_dfplot_bolt[j] > max_threshold)
				outliers_alpha_dfplot_bolt = alpha_dfplot_bolt[j].loc[condition_alpha_dfplot_bolt]
				all_outliers = set(outliers_alpha_dfplot.index.tolist() + outliers_alpha_dfplot_bolt.index.tolist())
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_bolt[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_bolt[j].loc[all_outliers], color = '#CA6627', s = 100, marker = '*')				
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_bolt[j], xerr = alpha_se_dfplot[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_bolt[j], yerr = alpha_se_dfplot_bolt[j], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$ BOLT-LMM, no PCs')
				ax[i].set_ylabel(r'$\hat{\alpha}$ BOLT-LMM with PCs')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if alpha_dfplot[j].loc[outlier] > alpha_dfplot_bolt[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_bolt[j].loc[outlier],label_dict[outlier],color = '#CA6627',horizontalalignment='right',verticalalignment='top')
					if alpha_dfplot[j].loc[outlier] < alpha_dfplot_bolt[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_bolt[j].loc[outlier],label_dict[outlier],color = '#CA6627')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/comparison_scatters/bolt.v.bolt.' + label + '.alpha.scatter.pdf')
			plt.clf()

	def plot_ascertainment_v_validation(self):
		for label in ['1kg.all', '1kg.eur']:
			fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (16,4))
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.half.1.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.half.1.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_df_val = pd.read_csv('../cache/alpha_matrices/plink.half.2.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df_val = pd.read_csv('../cache/alpha_matrices/plink.half.2.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot_val = alpha_df_val.loc[alpha_dfplot.index.tolist()]
			alpha_se_dfplot_val = alpha_se_df_val.loc[alpha_dfplot.index.tolist()]
			for i,j in zip([0,1,2,3],alpha_dfplot.columns):
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot[j].values)
				condition_alpha_dfplot = (alpha_dfplot[j] > max_threshold)
				outliers_alpha_dfplot = alpha_dfplot[j].loc[condition_alpha_dfplot]
				mu, sigma, min_threshold, max_threshold = self.estimate_gaussian_outliers(alpha_dfplot_val[j].values)
				condition_alpha_dfplot_val = (alpha_dfplot_val[j] > max_threshold)
				outliers_alpha_dfplot_val = alpha_dfplot_val[j].loc[condition_alpha_dfplot_val]
				all_outliers = set(outliers_alpha_dfplot.index.tolist() + outliers_alpha_dfplot_val.index.tolist())
				ax[i].title.set_text('p < ' + str(j))
				ax[i].scatter(alpha_dfplot[j].loc[~alpha_dfplot.index.isin(all_outliers)],alpha_dfplot_val[j].loc[~alpha_dfplot.index.isin(all_outliers)], color = '#CA6627', s = 35)
				ax[i].scatter(alpha_dfplot[j].loc[all_outliers],alpha_dfplot_val[j].loc[all_outliers], color = '#CA6627', s = 100, marker = '*')				
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_val[j], xerr = alpha_se_dfplot[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot_val[j], yerr = alpha_se_dfplot_val[j], color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].set_xlabel(r'$\hat{\alpha}$ ascertainment')
				ax[i].set_ylabel(r'$\hat{\alpha}$ validation')
				xseq = ax[i].get_xlim()
				ax[i].plot(xseq,xseq,color = 'grey',linestyle='--',alpha=0.75)
				xaxis_range = xseq[1]-xseq[0]
				for outlier in all_outliers:
					if alpha_dfplot[j].loc[outlier] > alpha_dfplot_val[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]-(xaxis_range*0.02),alpha_dfplot_val[j].loc[outlier],label_dict[outlier],color = '#CA6627',horizontalalignment='right',verticalalignment='top')
					if alpha_dfplot[j].loc[outlier] < alpha_dfplot_val[j].loc[outlier]:
						ax[i].text(alpha_dfplot[j].loc[outlier]+(xaxis_range*0.02),alpha_dfplot_val[j].loc[outlier],label_dict[outlier],color = '#CA6627')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/comparison_scatters/ascertainment.v.validation.' + label + '.alpha.scatter.pdf')
			plt.clf()



