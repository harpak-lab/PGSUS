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
import matplotlib.patches as mpatches
warnings.filterwarnings('ignore')

class thresholding_ascertainment(object):

	def __init__(self, label_dict, traits):
		self.threshold_palette = {'1.0':'#1b9e77', '0.001':'#7570b3', '1e-05':'#e7298a', '1e-08':'#d95f02'}
		self.shape_palette = {'1.0':'o', '0.001':'s', '1e-05':'^', '1e-08':'X'}
		self.spread = {'1.0':0.04, '0.001':0.02, '1e-05':-0.02, '1e-08':-0.04}
		self.label_dict = label_dict
		self.traits = traits

	def run(self):
		self.plot_alphas()
		
	def plot_alphas(self):
		best_best_all = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.best.clump.best.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		best_random_all = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.best.clump.random.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_best_all = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.random.clump.best.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_random_all = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.random.clump.random.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')

		best_best_all_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.best.clump.best.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		best_random_all_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.best.clump.random.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_best_all_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.random.clump.best.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_random_all_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.all.random.clump.random.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')

		best_best_eur = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.best.clump.best.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		best_random_eur = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.best.clump.random.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_best_eur = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.random.clump.best.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_random_eur = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.random.clump.random.snp.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')

		best_best_eur_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.best.clump.best.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		best_random_eur_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.best.clump.random.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_best_eur_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.random.clump.best.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		random_random_eur_se = pd.read_csv('../cache/thresholding_ascertainment_matrices/plink.wc.1kg.eur.random.clump.random.snp.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')

		for trait in self.traits:
			fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10,5), constrained_layout=True, sharey = True)
			best_best_all_plot = pd.DataFrame(np.vstack((np.array(best_best_all.loc[trait]),np.array(best_best_all_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			best_best_all_plot['color'] = best_best_all_plot.index.map(self.threshold_palette)
			best_best_all_plot['shape'] = best_best_all_plot.index.map(self.shape_palette)
			best_best_all_plot['x'] = 0.5 - best_best_all_plot.index.map(self.spread)

			best_random_all_plot = pd.DataFrame(np.vstack((np.array(best_random_all.loc[trait]),np.array(best_random_all_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			best_random_all_plot['color'] = best_random_all_plot.index.map(self.threshold_palette)
			best_random_all_plot['shape'] = best_random_all_plot.index.map(self.shape_palette)
			best_random_all_plot['x'] = 1 - best_random_all_plot.index.map(self.spread)

			random_best_all_plot = pd.DataFrame(np.vstack((np.array(random_best_all.loc[trait]),np.array(random_best_all_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			random_best_all_plot['color'] = random_best_all_plot.index.map(self.threshold_palette)
			random_best_all_plot['shape'] = random_best_all_plot.index.map(self.shape_palette)
			random_best_all_plot['x'] = 1.5 - random_best_all_plot.index.map(self.spread)

			random_random_all_plot = pd.DataFrame(np.vstack((np.array(random_random_all.loc[trait]),np.array(random_random_all_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			random_random_all_plot['color'] = random_random_all_plot.index.map(self.threshold_palette)
			random_random_all_plot['shape'] = random_random_all_plot.index.map(self.shape_palette)
			random_random_all_plot['x'] = 2 - random_random_all_plot.index.map(self.spread)

			for marker, d in best_best_all_plot.groupby('shape'):
				ax[0].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[0].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in best_random_all_plot.groupby('shape'):
				ax[0].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[0].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in random_best_all_plot.groupby('shape'):
				ax[0].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[0].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in random_random_all_plot.groupby('shape'):
				ax[0].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[0].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')

			best_best_eur_plot = pd.DataFrame(np.vstack((np.array(best_best_eur.loc[trait]),np.array(best_best_eur_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			best_best_eur_plot['color'] = best_best_eur_plot.index.map(self.threshold_palette)
			best_best_eur_plot['shape'] = best_best_eur_plot.index.map(self.shape_palette)
			best_best_eur_plot['x'] = 0.5 - best_best_eur_plot.index.map(self.spread)

			best_random_eur_plot = pd.DataFrame(np.vstack((np.array(best_random_eur.loc[trait]),np.array(best_random_eur_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			best_random_eur_plot['color'] = best_random_eur_plot.index.map(self.threshold_palette)
			best_random_eur_plot['shape'] = best_random_eur_plot.index.map(self.shape_palette)
			best_random_eur_plot['x'] = 1 - best_random_eur_plot.index.map(self.spread)

			random_best_eur_plot = pd.DataFrame(np.vstack((np.array(random_best_eur.loc[trait]),np.array(random_best_eur_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			random_best_eur_plot['color'] = random_best_eur_plot.index.map(self.threshold_palette)
			random_best_eur_plot['shape'] = random_best_eur_plot.index.map(self.shape_palette)
			random_best_eur_plot['x'] = 1.5 - random_best_eur_plot.index.map(self.spread)

			random_random_eur_plot = pd.DataFrame(np.vstack((np.array(random_random_eur.loc[trait]),np.array(random_random_eur_se.loc[trait]))), columns = ['1e-08','1e-05','0.001','1.0'], index = ['alpha','se']).T
			random_random_eur_plot['color'] = random_random_eur_plot.index.map(self.threshold_palette)
			random_random_eur_plot['shape'] = random_random_eur_plot.index.map(self.shape_palette)
			random_random_eur_plot['x'] = 2 - random_random_eur_plot.index.map(self.spread)

			for marker, d in best_best_eur_plot.groupby('shape'):
				ax[1].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[1].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in best_random_eur_plot.groupby('shape'):
				ax[1].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[1].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in random_best_eur_plot.groupby('shape'):
				ax[1].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = '')
				ax[1].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')
			for marker, d in random_random_eur_plot.groupby('shape'):
				ax[1].scatter(d['x'], d['alpha'], marker = marker,color = d.iloc[0]['color'], s = 10, label = d.index.tolist()[0])
				ax[1].errorbar(d['x'], d['alpha'], yerr = d['se'],marker = marker, color = d.iloc[0]['color'], label = '')

			ax[0].set_xticks([0.5,1,1.5,2])
			ax[0].set_xticklabels(['Best clump\nBest SNP','Best clump\nRandom SNP','Random clump\nBest SNP', 'Random clump\nRandom SNP'], fontsize = 14)
			ax[0].tick_params(labelrotation=45)
			
			ax[1].set_xticks([0.5,1,1.5,2])
			ax[1].set_xticklabels(['Best clump\nBest SNP','Best clump\nRandom SNP','Random clump\nBest SNP', 'Random clump\nRandom SNP'], fontsize = 14)
			ax[1].tick_params(labelrotation=45)

			one_patch = mpatches.Patch(color=self.threshold_palette['1.0'], label=r'$1$')
			three_patch = mpatches.Patch(color = self.threshold_palette['0.001'], label = r'$10^{-3}$')
			five_patch = mpatches.Patch(color = self.threshold_palette['1e-05'], label = r'$10^{-5}$')
			eight_patch = mpatches.Patch(color = self.threshold_palette['1e-08'], label = r'$10^{-8}$')
			ax[1].legend(handles=[one_patch, three_patch, five_patch,eight_patch], fontsize = 12)


			ax[0].axhline(y=1,xmin = 0, xmax = 2.25, linestyle='--', color = 'grey', alpha = 0.75)
			ax[1].axhline(y=1,xmin = 0, xmax = 2.25, linestyle='--', color = 'grey', alpha = 0.75)
			# fig.suptitle(self.label_dict[trait], size = 'large')
			fig.subplots_adjust(top=0.9)
			ax[0].set_title(self.label_dict[trait] + ', 1KG', size = 'medium')
			ax[1].set_title(self.label_dict[trait] + ', 1KG Europeans', size = 'medium')
			ax[0].set_ylabel(r'$\hat{\alpha}$', fontsize = 14)
			ax[1].set_ylabel(r'$\hat{\alpha}$', fontsize = 14)
			plt.tight_layout()
			sns.despine()
			plt.savefig('../figures/thresholding_ascertainment_plots/' + trait + '.pdf')
			plt.clf()









