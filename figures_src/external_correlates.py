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
from scipy.stats import spearmanr
import matplotlib as mpl
from matplotlib.lines import Line2D
warnings.filterwarnings('ignore')

class external_correlates(object):

	def __init__(self, analyses, label_dict):
		self.analyses = analyses
		self.label_dict = label_dict

	def run(self):
		if 'pc1' in self.analyses:
			self.plot_fvc()
		if 'pc_loadings' in self.analyses:
			self.plot_pcloadings()
		if 'ses_h2' in self.analyses:
			self.plot_ses_h2()
		if 'perm_v_normal_ses' in self.analyses:
			self.plot_perm_v_normal_ses()
		if 'component_scatters' in self.analyses:
			self.plot_component_scatters()

	def plot_perm_v_normal_ses(self):
		sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
			'hand_grip_strength','height','hip_circ','household_income','neuroticism_score',
			'overall_health','pack_years_smoking','pulse_rate','skin_color','waist_circ','years_schooling']
		for trait in sps_traits:
			df = pd.read_csv('../cache/misc_supplemental_data/se_comparison/' + trait + '.se.comparison.txt',sep = '\t')
			fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5,5))
			ax.scatter(np.array(df['normal_se']),np.array(df['EMP_SE']), alpha = 0.5, color = 'grey')
			ax.set_xlabel('Normal theory S.E.')
			ax.set_ylabel('Permutation S.E.')
			xseq = ax.get_xlim()
			ax.plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)
			ax.set_title(self.label_dict[trait])
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/external_correlates/se.comparison.' + trait + '.pdf')
			plt.clf()
	
	def plot_ses_h2(self):
		all_alphas = pd.read_csv('../cache/alpha_matrices/plink.wc.1kg.all.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep='\t').set_index('Unnamed: 0')
		eur_alphas = pd.read_csv('../cache/alpha_matrices/plink.wc.1kg.eur.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep='\t').set_index('Unnamed: 0')
		
		h2 = pd.read_csv('../cache/misc_supplemental_data/combined.trait.h2s.txt', sep = '\t').set_index('Unnamed: 0')
		ses = pd.read_csv('../cache/misc_supplemental_data/combined.trait.townsend.pearson.r.txt', sep = '\t').set_index('Unnamed: 0')
		fig, ax = plt.subplots(nrows = 1, ncols = 4, figsize = (20,5))
		thresholds = ['1e-08','1e-05','0.001','1.0']
		thresh_labels = [r'$1\times10^{-8}$',r'$1\times10^{-5}$',r'$1\times10^{-3}$','1']
		idx = np.arange(4)
		
		#make the h2 figure
		for thresh,thresh_label,i in zip(thresholds,thresh_labels,idx):
			ax[i].scatter(h2['ldsc_h2'],eur_alphas[thresh],color = '#CA6627')
			ax[i].set_xlabel(r'LDSC $h^2$',fontsize = 16)
			ax[i].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG Europeans',fontsize = 16)
			pearson_all = stats.pearsonr(h2['ldsc_h2'],eur_alphas[thresh])
			ax[i].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

		plt.show()
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/h2.alpha.correlations.pdf')
		plt.clf()

		fig, ax = plt.subplots(nrows = 2, ncols = 4, figsize = (20,10))
		#make the ses figure
		for thresh,i in zip(thresholds,idx):
			ax[0,i].scatter(ses['pearson_r'],all_alphas[thresh],color = '#CA6627')
			ax[0,i].set_xlabel('Trait correlation with\nTownsend deprivation index',fontsize = 16)
			ax[0,i].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG',fontsize = 16)
			pearson_all = stats.pearsonr(ses['pearson_r'],all_alphas[thresh])
			ax[0,i].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

			ax[1,i].scatter(ses['pearson_r'],eur_alphas[thresh],color = '#CA6627')
			ax[1,i].set_xlabel('Trait correlation with\nTownsend deprivation index',fontsize = 16)
			ax[1,i].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG Europeans',fontsize = 16)
			pearson_all = stats.pearsonr(ses['pearson_r'],eur_alphas[thresh])
			ax[1,i].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

		plt.show()
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/townsend.alpha.correlations.pdf')
		plt.clf()


	def plot_pcloadings(self):
		all_ukb = np.abs(pd.read_csv('../cache/misc_supplemental_data/loading.corr.all1kg.ukb.txt',sep = '\t').set_index('Unnamed: 0'))
		eur_ukb = np.abs(pd.read_csv('../cache/misc_supplemental_data/loading.corr.eur1kg.ukb.txt',sep = '\t').set_index('Unnamed: 0'))
		all_eur = np.abs(pd.read_csv('../cache/misc_supplemental_data/loading.corr.all1kg.eur1kg.txt',sep = '\t').set_index('Unnamed: 0'))

		fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (15,5))
		sns.heatmap(all_ukb, ax = ax[0], cmap='Blues', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[0].set_xlabel('UK Biobank',fontsize = 16)
		ax[0].set_xticks([i-0.5 for i in range(1,21)])
		ax[0].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=30,ha='center',fontsize = 7)
		ax[0].set_ylabel('1KG',fontsize = 16)
		ax[0].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 7)

		sns.heatmap(eur_ukb, ax = ax[1], cmap='Greens', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[1].set_xlabel('UK Biobank',fontsize = 16)
		ax[1].set_xticks([i-0.5 for i in range(1,21)])
		ax[1].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=30,ha='center',fontsize = 7)
		ax[1].set_ylabel('1KG Europeans',fontsize = 16)
		ax[1].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 7)

		sns.heatmap(all_eur, ax = ax[2], cmap='Purples', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[2].set_xlabel('1KG Europeans',fontsize = 16)
		ax[2].set_xticks([i-0.5 for i in range(1,21)])
		ax[2].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=30,ha='center',fontsize = 7)
		ax[2].set_ylabel('1KG',fontsize = 16)
		ax[2].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 7)

		plt.show()
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/pcloading.correlations.pdf')
		plt.clf()

	def plot_fvc(self):		
		pc_df = pd.read_csv('../cache/misc_supplemental_data/wc.first.pc.coordinates.txt',sep = '\t')[['FID','PC1']]
		pack_years_smoking_df = pd.read_csv('../cache/misc_supplemental_data/pheno_pack_years_smoking.txt',sep = '\t')[['FID','pack_years_smoking']]
		years_schooling_df = pd.read_csv('../cache/misc_supplemental_data/pheno_years_schooling.txt',sep = '\t')[['FID','years_schooling']]
		ids = pd.read_csv('../cache/misc_supplemental_data/wc.chr22.nosibs.nowithdrawals.sample', delim_whitespace = True).iloc[1:]
		ids = ids[['ID_1']].rename(columns={'ID_1':'FID'})

		smoking_plot_df = pack_years_smoking_df.merge(pc_df, on = 'FID', how = 'inner').merge(ids, on = 'FID', how = 'inner')
		schooling_plot_df = years_schooling_df.merge(pc_df, on = 'FID', how = 'inner').merge(ids, on = 'FID', how = 'inner')
		labels = [1,2,3,4,5,6,7,8,9,10]
		smoking_plot_df = smoking_plot_df.sort_values(by='PC1')
		tenth = smoking_plot_df.shape[0]/10
		smoking_plot = []
		for i in range(1,11):
			lower = (i-1)*tenth
			upper = ((i)*tenth)-1
			temp = smoking_plot_df.iloc[lower:upper]['pack_years_smoking'].mean()
			temp2 = smoking_plot_df.iloc[lower:upper]['PC1'].mean()
			smoking_plot.append([temp2,temp])

		schooling_plot_df = schooling_plot_df.sort_values(by='PC1').reset_index(drop=True)
		tenth = schooling_plot_df.shape[0]/10
		school_plot = []
		for i in range(1,11):
			lower = (i-1)*tenth
			upper = ((i)*tenth)-1
			temp = schooling_plot_df.iloc[lower:upper]['years_schooling'].mean()
			temp2 = schooling_plot_df.iloc[lower:upper]['PC1'].mean()
			school_plot.append([temp2,temp])

		smoking_plot = np.array(smoking_plot)
		schooling_plot = np.array(school_plot)
		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10,5))
		#first column where p<1 is always the x-axis
		ax[0].scatter(smoking_plot[:,0],smoking_plot[:,1], color = '#377eb8', s = 45)
		ax[0].set_xlabel(r'PC1')
		ax[0].set_ylabel(r'Pack Years of Smoking')
		smoking_fit = np.polyfit(smoking_plot_df['PC1'],smoking_plot_df['pack_years_smoking'],1)
		xmax = ax[0].get_xlim()[1]
		xmin = ax[0].get_xlim()[0]
		ax[0].plot([xmin,xmax],[xmin*smoking_fit[0] + smoking_fit[1],xmax*smoking_fit[0] + smoking_fit[1]])
		print('Pack Years of Smoking')
		print(pearsonr(smoking_plot_df['PC1'],smoking_plot_df['pack_years_smoking']))
		print(smoking_plot_df.shape[0])

		ax[1].scatter(schooling_plot[:,0],schooling_plot[:,1], color = '#4daf4a', s = 45)
		ax[1].set_xlabel(r'PC1')
		ax[1].set_ylabel(r'Years of Schooling')
		schooling_fit = np.polyfit(schooling_plot_df['PC1'],schooling_plot_df['years_schooling'],1)
		xmax = ax[0].get_xlim()[1]
		xmin = ax[0].get_xlim()[0]
		ax[1].plot([xmin,xmax],[xmin*schooling_fit[0] + schooling_fit[1],xmax*schooling_fit[0] + schooling_fit[1]], color = '#4daf4a')
		print('Years of Schooling')
		print(spearmanr(schooling_plot_df['PC1'],schooling_plot_df['years_schooling']))
		print(schooling_plot_df.shape[0])

		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/pc1.v.pheno.smoking.schooling.png',dpi=300)
		plt.clf()

	def plot_component_scatters(self):
		label = '1kg.all'
		alpha_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep = '\t').set_index('Unnamed: 0')
		alpha_se_df = pd.read_csv('../cache/alpha_matrices/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.alpha.se.mat.txt',sep = '\t').set_index('Unnamed: 0')
		threshold_palette_nondirect = {'1.0':'#e0905e', '0.001':'#d87333', '1e-05':'#b55b23', '1e-08':'#8a461b'}
		threshold_palette_sad = {'1.0':'#D97567', '0.001':'#CF4D3C', '1e-05':'#B53423', '1e-08':'#8D281B'}
		threshold_shapes = {'1.0':'s', '0.001':'p', '1e-05':'X', '1e-08':'*'}
		nondirect_df = pd.read_csv('../cache/alpha_matrices_nondirect/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.nondirect.variance.prop.mat.txt',sep = '\t').set_index('Unnamed: 0')
		sad_df = pd.read_csv('../cache/alpha_matrices_nondirect/plink.wc.' + label + '.sps23.aperm.1K.to.1M.v2.sad.variance.prop.mat.txt',sep = '\t').set_index('Unnamed: 0')
		alpha_dfplot = alpha_df.astype(float).sort_values(by='1e-05')
		alpha_se_df = alpha_se_df.loc[alpha_dfplot.index.tolist()]

		sad_df = sad_df.reset_index().melt(id_vars = ['Unnamed: 0'], value_vars = threshold_palette_sad.keys()).rename(columns = {'Unnamed: 0': 'traits', 'variable':'threshold', 'value':'sad'})
		nondirect_df = nondirect_df.reset_index().melt(id_vars = ['Unnamed: 0'], value_vars = threshold_palette_sad.keys()).rename(columns = {'Unnamed: 0': 'traits', 'variable':'threshold', 'value':'nondirect'})
		alpha_dfplot = alpha_dfplot.reset_index().melt(id_vars = ['Unnamed: 0'], value_vars = threshold_palette_sad.keys()).rename(columns = {'Unnamed: 0': 'traits', 'variable':'threshold', 'value':'alpha'})

		merged = pd.merge(left=sad_df, right=nondirect_df, how='left', left_on=['traits', 'threshold'], right_on=['traits', 'threshold'])
		merged = pd.merge(left=merged, right=alpha_dfplot, how='left', left_on=['traits', 'threshold'], right_on=['traits', 'threshold'])
		merged['sad_color'] = merged['threshold'].map(threshold_palette_sad)
		merged['nondirect_color'] = merged['threshold'].map(threshold_palette_nondirect)
		merged['shape'] = merged['threshold'].map(threshold_shapes)

		fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10,5))	
		for i,row in merged.iterrows():	
  
			ax[0].scatter(row['sad'],row['alpha'],color = row['sad_color'], marker = row['shape'], s = 45)
			ax[0].set_axisbelow(True)
			ax[0].set_ylabel('Isotropic inflation factor')
			ax[0].set_xlabel('SAD variance divided by non-error variance')

			ax[1].scatter(row['nondirect'],row['alpha'],color = row['nondirect_color'], marker = row['shape'], s = 45)
			ax[1].set_axisbelow(True)
			ax[1].set_ylabel('Isotropic inflation factor')
			ax[1].set_xlabel('Non-direct variance divided by\nnon-error variance')

		mypoints = [Line2D([0], [0], marker=threshold_shapes[threshold], markersize=10, markeredgecolor = threshold_palette_sad[threshold], markerfacecolor=threshold_palette_sad[threshold], linestyle='') for threshold in ['1.0','0.001','1e-05','1e-08']]
		mylabels = [r'$p$-value < ' + str(threshold) for threshold in ['1.0','0.001','1e-05','1e-08']]
		ax[0].legend(mypoints, mylabels)		
		mypoints = [Line2D([0], [0], marker=threshold_shapes[threshold], markersize=10, markeredgecolor = threshold_palette_nondirect[threshold], markerfacecolor=threshold_palette_nondirect[threshold], linestyle='') for threshold in ['1.0','0.001','1e-05','1e-08']]
		mylabels = [r'$p$-value < ' + str(threshold) for threshold in ['1.0','0.001','1e-05','1e-08']]
		ax[1].legend(mypoints, mylabels)

		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/fig.' + label + '.nonerror.scatters.pdf')
		plt.clf()
