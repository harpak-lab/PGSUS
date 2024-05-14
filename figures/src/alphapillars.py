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

class alphapillars(object):

	def __init__(self, analyses, label_dict):
		
		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		if 'wc' in self.analyses:
			self.plot_wc()
		if 'nopcs' in self.analyses:
			self.plot_nopcs()
		if 'halves' in self.analyses:
			self.plot_halves()
		if 'bolt' in self.analyses:
			self.plot_bolt()
		if 'nosingletons' in self.analyses:
			self.plot_nosingletons()
		if 'shuffle' in self.analyses:
			self.plot_shuffle()

	def plot_wc(self):
		for label in ['1kg.all', '1kg.eur']:	
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
			nsnp = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.nsnp.mat.txt',sep = '\t').set_index('Unnamed: 0')
			fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
						gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})
			for i,j in zip([0,2,4,6],alpha_dfplot.columns):
				ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
				ax[i].title.set_text('p < ' + str(j))
				ax[i].set_axisbelow(True)
				ax[i].set_xlabel(r'$\hat{\alpha}$')
				ax[i].grid(True)
				if i == 0:
					ax[i].set_yticks(alpha_dfplot['ycoordinate'])
					ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
				else:
					ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
					ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
				ax[i+1].set_ylim(ax[i].get_ylim())
				ax[i+1].axis('off')
				for k in alpha_dfplot.index.tolist():
					ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/alpha_pillars/wc.' + label + '.alphas.pdf')
			plt.clf()

	def plot_nopcs(self):
		for label in ['1kg.all', '1kg.eur']:
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.nopcs.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
			fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
						gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})
			for i,j in zip([0,2,4,6],alpha_dfplot.columns):
				ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
				ax[i].title.set_text('p < ' + str(j))
				ax[i].set_axisbelow(True)
				ax[i].set_xlabel(r'$\hat{\alpha}$')
				ax[i].grid(True)
				if i == 0:
					ax[i].set_yticks(alpha_dfplot['ycoordinate'])
					ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
				else:
					ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
					ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
				ax[i+1].set_ylim(ax[i].get_ylim())
				ax[i+1].axis('off')
				for k in alpha_dfplot.index.tolist():
					ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/alpha_pillars/wc.nopcs.' + label + '.alphas.pdf')
			plt.clf()

	def plot_halves(self):
		for half in ['1','2']:
			for label in ['1kg.all', '1kg.eur']:
				alpha_df = pd.read_csv('../cache/alpha_matrices/plink.half.' + half + '.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.half.' + half + '.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
				alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
				alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
				fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
							gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})
				for i,j in zip([0,2,4,6],alpha_dfplot.columns):
					ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
					ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
					ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
					ax[i].title.set_text('p < ' + str(j))
					ax[i].set_axisbelow(True)
					ax[i].set_xlabel(r'$\hat{\alpha}$')
					ax[i].grid(True)
					if i == 0:
						ax[i].set_yticks(alpha_dfplot['ycoordinate'])
						ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
					else:
						ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
						ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
					ax[i+1].set_ylim(ax[i].get_ylim())
					ax[i+1].axis('off')
					for k in alpha_dfplot.index.tolist():
						ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')
				sns.despine()
				plt.tight_layout()
				plt.savefig('../figures/alpha_pillars/half.' + half + '.' + label + '.alphas.pdf')
				plt.clf()

	def plot_bolt(self):
		for grm in ['nopcs','wpcs']:
			for label in ['1kg.all', '1kg.eur']:
				alpha_df = pd.read_csv('../cache/alpha_matrices/bolt.' + grm + '.' + label + '.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_se_df = pd.read_csv('../cache/alpha_matrices/bolt.' + grm + '.' + label + '.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
				alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
				alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
				alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
				fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
							gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})
				for i,j in zip([0,2,4,6],alpha_dfplot.columns):
					ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
					ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
					ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
					ax[i].title.set_text('p < ' + str(j))
					ax[i].set_axisbelow(True)
					ax[i].set_xlabel(r'$\hat{\alpha}$')
					ax[i].grid(True)
					if i == 0:
						ax[i].set_yticks(alpha_dfplot['ycoordinate'])
						ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
					else:
						ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
						ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
					ax[i+1].set_ylim(ax[i].get_ylim())
					ax[i+1].axis('off')
					for k in alpha_dfplot.index.tolist():
						ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')
				sns.despine()
				plt.tight_layout()
				plt.savefig('../figures/alpha_pillars/bolt.' + grm + '.' + label +'.alphas.pdf')
				plt.clf()
	
	def plot_nosingletons(self):
		for label in ['1kg.all', '1kg.eur']:
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.nosingletons.sps23.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.nosingletons.sps23.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
			
			fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
				gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})

			for i,j in zip([0,2,4,6],alpha_dfplot.columns):
				
				ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
				ax[i].title.set_text('p < ' + str(j))
				ax[i].set_axisbelow(True)
				ax[i].set_xlabel(r'$\hat{\alpha}$')
				ax[i].grid(True)
				
				if i == 0:
					ax[i].set_yticks(alpha_dfplot['ycoordinate'])
					ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
				else:
					ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
					ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
				
				ax[i+1].set_ylim(ax[i].get_ylim())
				ax[i+1].axis('off')
				for k in alpha_dfplot.index.tolist():
					ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')

			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/alpha_pillars/wc.nosingletons.' + label + '.alphas.pdf')
			plt.clf()

	def plot_shuffle(self):
		for label in ['1kg.all', '1kg.eur']:	
			alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.shuffle.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
			alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]
			alpha_dfplot['ycoordinate'] = [i for i in range(alpha_dfplot.shape[0],0,-1)]
			fig, ax = plt.subplots(nrows = 1, ncols = 8, figsize = (20,5), \
						gridspec_kw = {'width_ratios' :[4,0.25,4,0.25,4,0.25,4,0.25]})
			for i,j in zip([0,2,4,6],alpha_dfplot.columns):
				ax[i].vlines(1, ymin = 0.5, ymax =18.5, linestyles='dashed', color = 'black')
				ax[i].errorbar(alpha_dfplot[j],alpha_dfplot['ycoordinate'], xerr = alpha_se_df[j],color = '#CA6627', linestyle = '', capsize = 3)
				ax[i].scatter(alpha_dfplot[j],alpha_dfplot['ycoordinate'],color = '#CA6627', s = 45)
				ax[i].title.set_text('p < ' + str(j))
				ax[i].set_axisbelow(True)
				ax[i].set_xlabel(r'$\hat{\alpha}$')
				ax[i].grid(True)
				if i == 0:
					ax[i].set_yticks(alpha_dfplot['ycoordinate'])
					ax[i].set_yticklabels([self.label_dict[trait] for trait in alpha_dfplot.index.tolist()])
				else:
					ax[i].set_yticks([x+1 for x in range(alpha_dfplot.shape[0])])
					ax[i].set_yticklabels(['' for x in range(alpha_dfplot.shape[0])])
				ax[i+1].set_ylim(ax[i].get_ylim())
				ax[i+1].axis('off')
				for k in alpha_dfplot.index.tolist():
					ax[i+1].text(-3,alpha_dfplot.loc[k,'ycoordinate'] ,str(nsnp.loc[k][j]) + ' SNPs', verticalalignment = 'center')
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/alpha_pillars/wc.shuffle.' + label + '.alphas.pdf')
			plt.clf()

