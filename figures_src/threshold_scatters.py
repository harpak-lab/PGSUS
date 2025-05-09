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
from scipy.stats import spearmanr
import matplotlib as mpl
warnings.filterwarnings('ignore')

class threshold_scatters(object):

	def __init__(self, analyses, label_dict):
		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		if 'wc' in self.analyses:
			self.plot_wc()
		if 'nopcs' in self.analyses:
			self.plot_nopcs()
		if 'bolt' in self.analyses:
			self.plot_bolt()
		if 'shuffle' in self.analyses:
			self.plot_shuffle()

	def plot_wc(self):
		for label in ['1kg.all', '1kg.eur']:	
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			fig, ax = plt.subplots(nrows = 3, ncols = 3, figsize = (12,12))
			#first column where p<1 is always the x-axis
			ax[0,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['0.001'], color = '#CA6627', s = 45)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], yerr = alpha_se_dfplot['0.001'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['1.0'],alpha_dfplot['0.001'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[0,0].title.set_text(title)
			ax[0,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$', fontsize = 18)
			ax[0,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-3}$', fontsize = 18)
			xseq = ax[0,0].get_xlim()
			ax[0,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[1,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['1.0'],alpha_dfplot['1e-05'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[1,0].title.set_text(title)
			ax[1,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$', fontsize = 18)
			ax[1,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-5}$', fontsize = 18)
			xseq = ax[1,0].get_xlim()
			ax[1,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['1.0'],alpha_dfplot['1e-08'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[2,0].title.set_text(title)
			ax[2,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$', fontsize = 18)
			ax[2,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-8}$', fontsize = 18)
			xseq = ax[2,0].get_xlim()
			ax[2,0].plot(xseq,xseq,color = 'black',linestyle='--')
			#second column where p<0.001 is always the x-axis
			ax[1,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['0.001'],alpha_dfplot['1e-05'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[1,1].title.set_text(title)
			ax[1,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 10^{-3}$', fontsize = 18)
			ax[1,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-5}$', fontsize = 18)
			xseq = ax[1,1].get_xlim()
			ax[1,1].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['0.001'],alpha_dfplot['1e-08'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[2,1].title.set_text(title)
			ax[2,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 10^{-3}$', fontsize = 18)
			ax[2,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-8}$', fontsize = 18)
			xseq = ax[2,1].get_xlim()
			ax[2,1].plot(xseq,xseq,color = 'black',linestyle='--')
			#third column where p<1e-5 is always the x-axis
			ax[2,2].scatter(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			corr_coeff = spearmanr(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'])
			title = r'$\rho$ = ' + str(corr_coeff[0].round(3)) + '\n' + r'$p$-value $=$ ' + str(corr_coeff[1].round(3))
			ax[2,2].title.set_text(title)
			ax[2,2].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 10^{-5}$', fontsize = 18)
			ax[2,2].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 10^{-8}$', fontsize = 18)
			xseq = ax[2,2].get_xlim()
			ax[2,2].plot(xseq,xseq,color = 'black',linestyle='--')
			#set upper diagonal plots to whitespace
			ax[0,1].set_visible(False)
			ax[0,2].set_visible(False)
			ax[1,2].set_visible(False)
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/threshold_scatters/wc.' + label + '.alpha.threshold.scatter.pdf')
			plt.clf()
	
	def plot_nopcs(self):
		for label in ['1kg.all', '1kg.eur']:	
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			fig, ax = plt.subplots(nrows = 3, ncols = 3, figsize = (12,12))
			#first column where p<1 is always the x-axis
			ax[0,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['0.001'], color = '#CA6627', s = 45)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], yerr = alpha_se_dfplot['0.001'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[0,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[0,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			xseq = ax[0,0].get_xlim()
			ax[0,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[1,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[1,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			xseq = ax[1,0].get_xlim()
			ax[1,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[2,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,0].get_xlim()
			ax[2,0].plot(xseq,xseq,color = 'black',linestyle='--')
			#second column where p<0.001 is always the x-axis
			ax[1,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			ax[1,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			xseq = ax[1,1].get_xlim()
			ax[1,1].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			ax[2,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,1].get_xlim()
			ax[2,1].plot(xseq,xseq,color = 'black',linestyle='--')
			#third column where p<1e-5 is always the x-axis
			ax[2,2].scatter(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,2].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			ax[2,2].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,2].get_xlim()
			ax[2,2].plot(xseq,xseq,color = 'black',linestyle='--')
			#set upper diagonal plots to whitespace
			ax[0,1].set_visible(False)
			ax[0,2].set_visible(False)
			ax[1,2].set_visible(False)
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/threshold_scatters/wc.nopcs.' + label + '.alpha.threshold.scatter.pdf')
			plt.clf()

	def plot_bolt(self):
		for grm in ['nopcs','wpcs']:
			for label in ['1kg.all', '1kg.eur']:
				alpha_df = pd.read_csv('../cache/alpha_matrices/bolt.' + grm + '.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_se_df = pd.read_csv('../cache/alpha_matrices/bolt.' + grm + '.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
				alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
				fig, ax = plt.subplots(nrows = 3, ncols = 3, figsize = (12,12))
				#first column where p<1 is always the x-axis
				ax[0,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['0.001'], color = '#CA6627', s = 45)
				ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], yerr = alpha_se_dfplot['0.001'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[0,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
				ax[0,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
				xseq = ax[0,0].get_xlim()
				ax[0,0].plot(xseq,xseq,color = 'black',linestyle='--')
				ax[1,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
				ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[1,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
				ax[1,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
				xseq = ax[1,0].get_xlim()
				ax[1,0].plot(xseq,xseq,color = 'black',linestyle='--')
				ax[2,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
				ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
				ax[2,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
				xseq = ax[2,0].get_xlim()
				ax[2,0].plot(xseq,xseq,color = 'black',linestyle='--')
				#second column where p<0.001 is always the x-axis
				ax[1,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
				ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[1,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
				ax[1,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
				xseq = ax[1,1].get_xlim()
				ax[1,1].plot(xseq,xseq,color = 'black',linestyle='--')
				ax[2,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
				ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
				ax[2,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
				xseq = ax[2,1].get_xlim()
				ax[2,1].plot(xseq,xseq,color = 'black',linestyle='--')
				#third column where p<1e-5 is always the x-axis
				ax[2,2].scatter(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
				ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
				ax[2,2].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
				ax[2,2].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
				xseq = ax[2,2].get_xlim()
				ax[2,2].plot(xseq,xseq,color = 'black',linestyle='--')
				#set upper diagonal plots to whitespace
				ax[0,1].set_visible(False)
				ax[0,2].set_visible(False)
				ax[1,2].set_visible(False)
				sns.despine()
				plt.tight_layout()
				plt.savefig('../figures/threshold_scatters/bolt.' + grm + '.' + label + '.alpha.threshold.scatter.pdf')
				plt.clf()

	def plot_shuffle(self):
		for label in ['1kg.all', '1kg.eur']:	
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_dfplot = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			fig, ax = plt.subplots(nrows = 3, ncols = 3, figsize = (12,12))
			#first column where p<1 is always the x-axis
			ax[0,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['0.001'], color = '#CA6627', s = 45)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[0,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['0.001'], yerr = alpha_se_dfplot['0.001'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[0,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[0,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			xseq = ax[0,0].get_xlim()
			ax[0,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[1,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[1,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			xseq = ax[1,0].get_xlim()
			ax[1,0].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,0].scatter(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['1.0'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,0].errorbar(alpha_dfplot['1.0'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,0].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1$')
			ax[2,0].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,0].get_xlim()
			ax[2,0].plot(xseq,xseq,color = 'black',linestyle='--')
			#second column where p<0.001 is always the x-axis
			ax[1,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], color = '#CA6627', s = 45)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-05'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[1,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			ax[1,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			xseq = ax[1,1].get_xlim()
			ax[1,1].plot(xseq,xseq,color = 'black',linestyle='--')
			ax[2,1].scatter(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,1].errorbar(alpha_dfplot['0.001'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-08'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,1].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 0.001$')
			ax[2,1].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,1].get_xlim()
			ax[2,1].plot(xseq,xseq,color = 'black',linestyle='--')
			#third column where p<1e-5 is always the x-axis
			ax[2,2].scatter(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], color = '#CA6627', s = 45)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], xerr = alpha_se_dfplot['0.001'],color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,2].errorbar(alpha_dfplot['1e-05'],alpha_dfplot['1e-08'], yerr = alpha_se_dfplot['1e-05'], color = '#CA6627', linestyle = '', capsize = 3)
			ax[2,2].set_xlabel(r'$\hat{\alpha}$ $p$-value $< 1e-05$')
			ax[2,2].set_ylabel(r'$\hat{\alpha}$ $p$-value $< 1e-08$')
			xseq = ax[2,2].get_xlim()
			ax[2,2].plot(xseq,xseq,color = 'black',linestyle='--')
			#set upper diagonal plots to whitespace
			ax[0,1].set_visible(False)
			ax[0,2].set_visible(False)
			ax[1,2].set_visible(False)
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/threshold_scatters/wc.shuffle.' + label + '.alpha.threshold.scatter.pdf')
			plt.clf()

