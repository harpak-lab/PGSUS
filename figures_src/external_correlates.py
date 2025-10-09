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
import statsmodels.formula.api as smf
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
		if 'sib_se_comparison' in self.analyses:
			self.sib_se_comparison()
		if 'component_scatters' in self.analyses:
			self.plot_component_scatters()
		if 'bjk_se' in self.analyses:
			self.bjk_se()

	def bjk_se(self):
		df = pd.read_csv('../cache/misc_supplemental_data/se_comparison/bjk.se.pheno.var.data.txt', sep = '\t').set_index('Unnamed: 0')
		fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (15,5))

		ax[0].scatter(np.sqrt(np.array(df['pop_pheno_var'])),np.sqrt(np.array(df['pop_sib_var'])), alpha = 0.5, color = 'grey')
		ax[0].set_xlabel('S.D. of population phenotypes')
		ax[0].set_ylabel('S.D. of sibling phenotypes')
		xseq = ax[0].get_xlim()
		x = np.sqrt(np.array(df['pop_pheno_var']))
		y = np.sqrt(np.array(df['pop_sib_var']))
		n = [self.label_dict[j] for j in df.index.tolist()]

		for i, txt in enumerate(n):
			ax[0].text(x[i], y[i], txt, horizontalalignment = 'center')

		ax[0].plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)

		ax[1].scatter(np.sqrt(np.array(df['pop_pheno_var'])),np.array(df['bjk_se ~ EMP_SE']), alpha = 0.5, color = 'grey')
		ax[1].set_xlabel('S.D. of population phenotypes')
		ax[1].set_ylabel('Slope of Block jackknife S.E.\n regressed on permutation S.E.')
		# xseq = ax[1].get_xlim()
		# ax[1].plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)

		ax[2].scatter(np.sqrt(np.array(df['pop_sib_var'])),np.array(df['bjk_se ~ EMP_SE']), alpha = 0.5, color = 'grey')
		ax[2].set_xlabel('S.D. of sibling phenotypes')
		ax[2].set_ylabel('Slope of Block jackknife S.E.\n regressed on permutation S.E.')

		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/bjk.comparisons.pdf')
		plt.clf()

	def sib_se_comparison(self):
		sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
			'hand_grip_strength','height','hip_circ','household_income','neuroticism_score',
			'overall_health','pack_years_smoking','pulse_rate','skin_color','waist_circ','years_schooling']
		sps_traits = ['height','bmi', 'years_schooling']
		for trait in sps_traits:
			df = pd.read_csv('../cache/misc_supplemental_data/plink.wc.'+trait+'.aperm.1K.to.1M.sib.preproc.txt',sep = '\t')
			normal_se = pd.read_csv('../cache/misc_supplemental_data/se_comparison/'+trait+'.se.comparison.txt',sep = '\t')
			df['normal_se'] = normal_se['normal_se']
			bjk = pd.read_csv('../cache/misc_supplemental_data/'+trait+'.bjk.se.sib.estimates.txt',sep = '\t')
			df = df.merge(bjk[['SNP']], on = 'SNP', how = 'inner')

			fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (15,5))
			ax[0].scatter(np.array(df['normal_se']),np.array(df['EMP_SE']), alpha = 0.5, color = 'grey')
			ax[0].set_xlabel('Normal theory S.E.')
			ax[0].set_ylabel('Permutation S.E.')
			xseq = ax[0].get_xlim()
			ax[0].plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)
			ax[0].set_title(self.label_dict[trait])
			
			ax[1].scatter(np.array(df['normal_se']),np.array(bjk['bjk_se']), alpha = 0.5, color = 'grey')
			ax[1].set_xlabel('Normal theory S.E.')
			ax[1].set_ylabel('Block jackknife S.E.')
			xseq = ax[1].get_xlim()
			ax[1].plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)
			ax[1].set_title(self.label_dict[trait])
			

			model = sm.OLS(bjk['bjk_se'],df['EMP_SE']).fit()
			print(trait,1./model.params[0])
			ax[2].scatter(np.array(df['EMP_SE']), np.array(bjk['bjk_se']),alpha = 0.5, color = 'grey')
			ax[2].set_ylabel('Block jackknife S.E.')
			ax[2].set_xlabel('Permutation S.E.')
			xseq = ax[2].get_xlim()
			ax[2].plot(xseq,xseq,color = 'black',linestyle='--',alpha=0.75)
			ax[2].set_title(self.label_dict[trait])
			
			sns.despine()
			plt.tight_layout()
			plt.savefig('../figures/external_correlates/se.comparison.' + trait + '.png')
			plt.clf()
	
	def plot_ses_h2(self):
		all_alphas = pd.read_csv('../cache/alpha_matrices/plink.wc.1kg.all.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep='\t').set_index('Unnamed: 0')
		eur_alphas = pd.read_csv('../cache/alpha_matrices/plink.wc.1kg.eur.sps23.aperm.1K.to.1M.v2.alpha.mat.txt',sep='\t').set_index('Unnamed: 0')
		
		h2 = pd.read_csv('../cache/misc_supplemental_data/combined.trait.h2s.txt', sep = '\t').set_index('Unnamed: 0')
		ses = pd.read_csv('../cache/misc_supplemental_data/combined.trait.townsend.pearson.r.txt', sep = '\t').set_index('Unnamed: 0')
		fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (10,10))
		thresholds = ['1e-08','1e-05','0.001','1.0']
		thresh_labels = [r'$10^{-8}$',r'$10^{-5}$',r'$10^{-3}$','1']
		idx = [[0,0],[0,1],[1,0],[1,1]]
		print(h2['ldsc_h2'])
		print(eur_alphas['1.0'])
		#make the h2 figure
		for thresh,thresh_label,i in zip(thresholds,thresh_labels,idx):
			ax[i[0],i[1]].scatter(h2['ldsc_h2'],eur_alphas[thresh],color = '#CA6627')
			ax[i[0],i[1]].set_xlabel(r'LDSC $h^2$',fontsize = 16)
			ax[i[0],i[1]].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG Europeans',fontsize = 16)
			pearson_all = stats.pearsonr(h2['ldsc_h2'],eur_alphas[thresh])
			ax[i[0],i[1]].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

		plt.show()
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/h2.alpha.correlations.pdf')
		plt.clf()

		fig, ax = plt.subplots(nrows = 4, ncols = 2, figsize = (10,20))
		#make the ses figure
		idx = [0,1,2,3]
		for thresh,thresh_label,i in zip(thresholds,thresh_labels,idx):
			ax[i,0].scatter(ses['pearson_r'],all_alphas[thresh],color = '#CA6627')
			ax[i,0].set_xlabel('Trait correlation with\nTownsend deprivation index',fontsize = 16)
			ax[i,0].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG',fontsize = 16)
			pearson_all = stats.pearsonr(ses['pearson_r'],all_alphas[thresh])
			ax[i,0].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

			ax[i,1].scatter(ses['pearson_r'],eur_alphas[thresh],color = '#CA6627')
			ax[i,1].set_xlabel('Trait correlation with\nTownsend deprivation index',fontsize = 16)
			ax[i,1].set_ylabel(r'$\hat{\alpha}$ at $p$-value $<$ ' + thresh_label + '\n1KG Europeans',fontsize = 16)
			pearson_all = stats.pearsonr(ses['pearson_r'],eur_alphas[thresh])
			ax[i,1].set_title(r'Pearson r $=$ ' + str(pearson_all[0].round(4)) + '\n' + r'$p$-value $=$ ' + str(pearson_all[1].round(4)) ,fontsize = 16)

		plt.show()
		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/townsend.alpha.correlations.pdf')
		plt.clf()

	def plot_pcloadings(self):
		ukb_loadings = pd.read_csv('../cache/1kg_pc_data/supp.fig.ukb.loadings',delim_whitespace = True)
		ukb_loadings.columns = ['SNP','A1'] + ['ukbPC' + str(i) for i in range(1,21)]
		ukb_loadings = ukb_loadings[['SNP'] + ['ukbPC' + str(i) for i in range(1,21)]]
		all_loadings = pd.read_csv('../cache/1kg_pc_data/supp.fig.all1000G.loadings',delim_whitespace = True)
		all_loadings.columns = ['SNP','A1'] + ['allPC' + str(i) for i in range(1,21)]
		all_loadings = all_loadings[['SNP'] + ['allPC' + str(i) for i in range(1,21)]]
		eur_loadings = pd.read_csv('../cache/1kg_pc_data/supp.fig.eur1000G.loadings',delim_whitespace = True)
		eur_loadings.columns = ['SNP','A1'] + ['eurPC' + str(i) for i in range(1,21)]
		eur_loadings = eur_loadings[['SNP'] + ['eurPC' + str(i) for i in range(1,21)]]


		all_ukb = ukb_loadings.merge(all_loadings, on = 'SNP', how = 'inner')
		all_ukb = np.abs(all_ukb[all_ukb.columns[1:]].corr().loc[['ukbPC' + str(i) for i in range(1,21)]][['allPC' + str(i) for i in range(1,21)]])
		eur_ukb = ukb_loadings.merge(eur_loadings, on = 'SNP', how = 'inner')
		eur_ukb = np.abs(eur_ukb[eur_ukb.columns[1:]].corr().loc[['ukbPC' + str(i) for i in range(1,21)]][['eurPC' + str(i) for i in range(1,21)]])
		all_eur = all_loadings.merge(eur_loadings, on = 'SNP', how = 'inner')
		all_eur = np.abs(all_eur[all_eur.columns[1:]].corr().loc[['allPC' + str(i) for i in range(1,21)]][['eurPC' + str(i) for i in range(1,21)]])

		# eur_ukb = np.abs(pd.read_csv('../cache/misc_supplemental_data/loading.corr.eur1kg.ukb.txt',sep = '\t').set_index('Unnamed: 0'))
		# all_eur = np.abs(pd.read_csv('../cache/misc_supplemental_data/loading.corr.all1kg.eur1kg.txt',sep = '\t').set_index('Unnamed: 0'))

		fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (12,10))
		sns.heatmap(all_ukb, ax = ax[0,0], cmap='Blues', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[0,0].set_xlabel('UK Biobank',fontsize = 16)
		ax[0,0].set_xticks([i-0.5 for i in range(1,21)])
		ax[0,0].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=90,ha='center',fontsize = 14)
		ax[0,0].set_ylabel('1KG',fontsize = 16)
		ax[0,0].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 14)

		sns.heatmap(eur_ukb, ax = ax[0,1], cmap='Greens', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[0,1].set_xlabel('UK Biobank',fontsize = 16)
		ax[0,1].set_xticks([i-0.5 for i in range(1,21)])
		ax[0,1].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=90,ha='center',fontsize = 14)
		ax[0,1].set_ylabel('1KG Europeans',fontsize = 16)
		ax[0,1].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 14)

		sns.heatmap(all_eur, ax = ax[1,0], cmap='Purples', fmt = '.2f', cbar_kws={'label': '|r| between loadings'})
		ax[1,0].set_xlabel('1KG Europeans',fontsize = 16)
		ax[1,0].set_xticks([i-0.5 for i in range(1,21)])
		ax[1,0].set_xticklabels(['PC' + str(i) for i in range(1,21)],rotation=90,ha='center',fontsize = 14)
		ax[1,0].set_ylabel('1KG',fontsize = 16)
		ax[1,0].set_yticklabels(['PC' + str(i) for i in range(1,21)], fontsize = 14)

		ax[1,1].set_visible(False)
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
			ax[0].set_ylabel('Isotropic inflation factor', fontsize = 12)
			ax[0].set_xlabel('SAD variance divided by non-error variance', fontsize = 12)

			ax[1].scatter(row['nondirect'],row['alpha'],color = row['nondirect_color'], marker = row['shape'], s = 45)
			ax[1].set_axisbelow(True)
			ax[1].set_ylabel('Isotropic inflation factor', fontsize = 12)
			ax[1].set_xlabel('Non-direct variance divided by\nnon-error variance', fontsize = 12)

		mypoints = [Line2D([0], [0], marker=threshold_shapes[threshold], markersize=10, markeredgecolor = threshold_palette_sad[threshold], markerfacecolor=threshold_palette_sad[threshold], linestyle='') for threshold in ['1.0','0.001','1e-05','1e-08']]
		mylabels = [r'$p$-value < ' + str(threshold) for threshold in ['1',r'$10^{-3}$',r'$10^{-5}$',r'$10^{-8}$']]
		ax[0].legend(mypoints, mylabels, fontsize = 12)		
		mypoints = [Line2D([0], [0], marker=threshold_shapes[threshold], markersize=10, markeredgecolor = threshold_palette_nondirect[threshold], markerfacecolor=threshold_palette_nondirect[threshold], linestyle='') for threshold in ['1.0','0.001','1e-05','1e-08']]
		mylabels = [r'$p$-value < ' + str(threshold) for threshold in ['1',r'$10^{-3}$',r'$10^{-5}$',r'$10^{-8}$']]
		ax[1].legend(mypoints, mylabels, fontsize = 12)

		sns.despine()
		plt.tight_layout()
		plt.savefig('../figures/external_correlates/fig.' + label + '.nonerror.scatters.pdf')
		plt.clf()
