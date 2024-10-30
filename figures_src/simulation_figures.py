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

class simulation_plots(object):

	def __init__(self, analyses, label_dict):
		self.analyses = analyses
		self.label_dict = label_dict
		self.sad_palette = {'sad':'#ca3a27', 'direct':'#4B9C79', 'covar':'#D1BA41'}

	def run(self):
		if 'null' in self.analyses:
			self.plot_null()
		if 'thresholds' in self.analyses:
			self.plot_thresholds()
		if 'all_loci' in self.analyses:
			self.plot_allloci()
		if 'component_error' in self.analyses:
			self.plot_component_error()
		
	def plot_null(self):
		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10,5))
		arch_dict={'neutral':'Neutral','maf':'MAF-Dependent'}
		for k,architecture in enumerate(['neutral','maf']):

			alpha_01 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.1.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_01 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.1.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_02 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.2.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_02 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.2.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_05 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.5.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_05 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.5.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_08 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.8.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_08 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.0.8.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_01 = alpha_01[~(alpha_01==0.).any(axis=1)]
			alpha_plot_01 = alpha_01.mean(axis=0)
			alpha_se_01 = alpha_01.std(axis=0)
			alpha_plot_02 = alpha_02.mean(axis=0)
			alpha_se_02 = alpha_02.std(axis=0)
			alpha_plot_05 = alpha_05.mean(axis=0)
			alpha_se_05 = alpha_05.std(axis=0)
			alpha_plot_08 = alpha_08.mean(axis=0)
			alpha_se_08 = alpha_08.std(axis=0)

			ax[k].axhline(1,-0.25,3.25, linestyle = '--',color = 'black', label = r'$\alpha=1$')
			ax[k].scatter([i-0.2 for i in range(len(alpha_plot_01.index.tolist()))], alpha_plot_01, marker = '.', color = '#e7298a', s = 100, label = '0.1')
			ax[k].errorbar([i-0.2 for i in range(len(alpha_plot_01.index.tolist()))], alpha_plot_01, yerr=alpha_se_01, marker = '', linestyle = '',color = '#e7298a')
			ax[k].scatter([i-0.1 for i in range(len(alpha_plot_02.index.tolist()))], alpha_plot_02, marker = '.', color = '#1b9e77', s = 100, label = '0.2')
			ax[k].errorbar([i-0.1 for i in range(len(alpha_plot_02.index.tolist()))], alpha_plot_02, yerr=alpha_se_02, marker = '', linestyle = '',color = '#1b9e77')
			ax[k].scatter([i for i in range(len(alpha_plot_05.index.tolist()))], alpha_plot_05, marker = '.', color = '#d95f02', s = 100,  label = '0.5')
			ax[k].errorbar([i for i in range(len(alpha_plot_05.index.tolist()))], alpha_plot_05, yerr=alpha_se_05, marker = '', linestyle = '',color = '#d95f02')
			ax[k].scatter([i+0.1 for i in range(len(alpha_plot_08.index.tolist()))], alpha_plot_08, marker = '.', color = '#7570b3', s = 100,  label = '0.8')
			ax[k].errorbar([i+0.1 for i in range(len(alpha_plot_08.index.tolist()))], alpha_plot_08, yerr=alpha_se_08, marker = '', linestyle = '',color = '#7570b3')
			
			ax[k].set_xticks([0,1,2],['0.25','0.5','0.75'])
			ax[k].set_xlabel('Environmental covariance between siblings')
			ax[k].title.set_text(arch_dict[architecture] + ' genetic architecture')
		
		plt.suptitle('Average estimates of the isotropic inflation\nunder null models of genetic architecture')
		plt.legend(loc='best')
		plt.tight_layout()
		sns.despine()
		plt.savefig('../figures/simulations/null.architectures.pdf')

	def plot_thresholds(self):
		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (12,6),sharey=True)
		arch_dict={'neutral':'Neutral','maf':'MAF-Dependent'}
		for k,architecture in enumerate(['neutral','maf']):

			alpha_1 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.1.0.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_1 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.1.0.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_3 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.001.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_3 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.001.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_5 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.00001.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_5 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.00001.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_8 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.00000001.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se_8 = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.0.00000001.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_plot_1 = alpha_1.mean(axis=0)
			alpha_plot_3 = alpha_3.mean(axis=0)
			alpha_plot_5 = alpha_5.mean(axis=0)
			alpha_plot_8 = alpha_8.mean(axis=0)
			
			alpha_se_plot_1 = alpha_1.sem(axis=0)
			alpha_se_plot_3 = alpha_3.sem(axis=0)
			alpha_se_plot_5 = alpha_5.sem(axis=0)
			alpha_se_plot_8 = alpha_8.sem(axis=0)

			ax[k].axhline(1,-0.25,3.25, linestyle = '--',color = 'black', label = r'$\alpha=1$')
			ax[k].plot([], [], ' ', label="Ascertainment\n" + r"$p$-value")
			ax[k].plot([i-0.2 for i in range(len(alpha_plot_1.index.tolist()))], alpha_plot_1, marker = '.', color = '#e7298a', label = '1')
			ax[k].errorbar([i-0.2 for i in range(len(alpha_plot_1.index.tolist()))], alpha_plot_1, yerr=alpha_se_plot_1, marker = '', linestyle = '',color = '#e7298a')
			ax[k].plot([i-0.1 for i in range(len(alpha_plot_3.index.tolist()))], alpha_plot_3, marker = '.', color = '#1b9e77', label = r'$10^{-3}$')
			ax[k].errorbar([i-0.1 for i in range(len(alpha_plot_3.index.tolist()))], alpha_plot_3, yerr=alpha_se_plot_3, marker = '', linestyle = '',color = '#1b9e77')
			ax[k].plot([i for i in range(len(alpha_plot_5.index.tolist()))], alpha_plot_5, marker = '.', color = '#d95f02',  label = r'$10^{-5}$')
			ax[k].errorbar([i for i in range(len(alpha_plot_5.index.tolist()))], alpha_plot_5, yerr=alpha_se_plot_5, marker = '', linestyle = '',color = '#d95f02')
			ax[k].plot([i+0.1 for i in range(len(alpha_plot_8.index.tolist()))], alpha_plot_8, marker = '.', color = '#7570b3',  label = r'$10^{-8}$')
			ax[k].errorbar([i+0.1 for i in range(len(alpha_plot_8.index.tolist()))], alpha_plot_8, yerr=alpha_se_plot_8, marker = '', linestyle = '',color = '#7570b3')
			
			ax[k].set_xticks([0,1,2,3])
			ax[k].set_xticklabels(['0.1','0.2','0.5','0.8'])
			ax[k].set_xlabel('Environmental variance')
			ax[k].title.set_text(arch_dict[architecture] + ' genetic architecture')

		plt.suptitle('Average estimates of the isotropic inflation\nusing a range of ascertainment threholds')
		plt.legend(bbox_to_anchor = (1.03,0.65))
		plt.tight_layout()
		sns.despine()
		fig.subplots_adjust(top=0.85)
		plt.savefig('../figures/simulations/null.architectures.ascertainment.thresholds.pdf')


	def plot_allloci(self):
		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (12,6),sharey=True)
		arch_dict={'neutral':'Neutral','maf':'MAF-Dependent'}
		for k,architecture in enumerate(['neutral','maf']):

			alpha = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.0.8.nsib.20000.nstandard.50000.ntraitloci.400.alpha.txt',sep = '\t').set_index('Unnamed: 0')
			alpha_se = pd.read_csv('../cache/simulation_output_cache/null.' + architecture + '.architecture.env.var.0.8.nsib.20000.nstandard.50000.ntraitloci.400.alpha.se.txt',sep = '\t').set_index('Unnamed: 0')

			alpha_plot = alpha.mean(axis=0)
			alpha_se_plot = alpha.sem(axis=0)

			ax[k].axhline(1,-0.25,4.25, linestyle = '--',color = 'black', label = r'$\alpha=1$')
			ax[k].plot([i for i in range(len(alpha_plot.index.tolist()))], alpha_plot, marker = '.', color = '#CA6627',)
			ax[k].errorbar([i for i in range(len(alpha_plot.index.tolist()))], alpha_plot, yerr=alpha_se_plot, marker = '', linestyle = '',color = '#CA6627',)
			ax[k].set_xticks([0,1,2,3,4])
			ax[k].set_xticklabels(['0.01','0.1','0.2','0.5','0.8'])
			ax[k].set_xlabel('Environmental variance')
			ax[k].title.set_text(arch_dict[architecture] + ' genetic architecture')

		plt.suptitle('Average estimates of the isotropic inflation\nusing a range of environmental variance and 100% causal loci')
		plt.tight_layout()
		sns.despine()
		fig.subplots_adjust(top=0.85)
		plt.savefig('../figures/simulations/null.architectures.all.loci.pdf')


	def plot_component_error(self):
		env_vars = ['0.01','0.1','0.2','0.5','0.8']
		for architecture in ['neutral','maf']:
			mean_alpha_df = pd.read_csv('../cache/simulation_output_cache/mean.alpha.architecture.' + architecture + '.400.loci.txt',sep = '\t').set_index('Unnamed: 0')
			print(mean_alpha_df.index.tolist())

			for env_var in env_vars:
				df = pd.read_csv('../cache/simulation_output_cache/decomp.error.architecture.' + architecture + '.env.var.' + env_var + '.400.loci.txt',sep = '\t').set_index('Unnamed: 0')
				standard_error = df.loc['standard']
				sibling_error = df.loc['sib']

				fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (10,5))
		
				ax[0].scatter([j+101 for j,k in enumerate(standard_error)],standard_error,alpha = 0.5, color = 'grey')
				ax[0].set_ylabel(r"$c_{D i}'-\hat{c}_{D i}'$ (Standard)")
				ax[0].set_xlabel("Principal component")
				ax[0].axhline(y=0, color = 'black', linestyle = '--')
				
				ax[1].scatter([j+101 for j,k in enumerate(sibling_error)],sibling_error,alpha = 0.5, color = 'grey')
				ax[1].set_ylabel(r"$c_{D i}-\hat{c}_{D i}$ (Sibling)")
				ax[1].set_xlabel("Principal component")
				ax[1].axhline(y=0, color = 'black', linestyle = '--')
				print(mean_alpha_df)
				fig.suptitle(r'Average $\hat{\alpha}=$ ' + str(mean_alpha_df.loc[float(env_var)][0].round(6)))
				sns.despine()
				plt.tight_layout()
				plt.savefig('../figures/simulations/decomp.scatter.error.architecture.' + architecture + '.env.var.' + env_var + '.400.loci.pdf')






