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
import xlwt
warnings.filterwarnings('ignore')

class trait_components(object):

	def __init__(self, analyses, label_dict, traits):

		self.analyses = analyses
		self.label_dict = label_dict
		self.sps_traits = traits
		self.thresholds = [thresh for thresh in ['1.0','0.001','1e-05','1e-08']]

		self.summarydf_sad = pd.DataFrame(np.zeros((len(self.sps_traits), len(self.thresholds))), index = self.sps_traits, columns = self.thresholds)
		self.summarydf_direct = pd.DataFrame(np.zeros((len(self.sps_traits), len(self.thresholds))), index = self.sps_traits, columns = self.thresholds)
		self.summarydf_covar = pd.DataFrame(np.zeros((len(self.sps_traits), len(self.thresholds))), index = self.sps_traits, columns = self.thresholds)

	def run(self):
		if 'wc' in self.analyses:
			self.plot_wc()
		if 'nopcs' in self.analyses:
			self.plot_nopcs()
		if 'bolt' in self.analyses:
			self.plot_bolt()
		if 'maf01' in self.analyses:
			self.plot_maf01()
		if '1kg.pcs.only' in self.analyses:
			self.plot_1kg_pcs_only()
		if 'ukb.and.1kg.pcs' in self.analyses:
			self.plot_ukb_and_1kg()

	def plot_summaries(self,label,analysis):
		direct_colors = {'1.0':'#3c7c60','0.001':'#4b9c79','1e-05':'#6eaf93','1e-08':'#93c3ae'}
		sad_colors = {'1.0':'#a12e1f','0.001':'#ca3a27','1e-05':'#d46152','1e-08':'#df887d'}
		covar_colors = {'1.0':'#a79434','0.001':'#d1ba41','1e-05':'#dac766','1e-08':'#e3d58d'}
		legend_palette = {'1.0':'#666666','0.001':'#808080','1e-05':'#999999','1e-08':'#b2b2b2'}
		fig,ax = plt.subplots(1,3,figsize = (10,9))

		direct_plot = pd.melt(self.summarydf_direct.reset_index(), id_vars=['index'])
		covar_plot = pd.melt(self.summarydf_covar.reset_index(), id_vars=['index'])
		sad_plot = pd.melt(self.summarydf_sad.reset_index(), id_vars=['index'])

		sns.barplot(data = direct_plot, x ='value', y='index', hue = 'variable', orient = 'h', ax = ax[0], hue_order = ['1.0','0.001','1e-05','1e-08'], palette=direct_colors)
		sns.barplot(data = covar_plot, x ='value', y='index', hue = 'variable', orient = 'h', ax = ax[1], hue_order = ['1.0','0.001','1e-05','1e-08'], palette=covar_colors)
		sns.barplot(data = sad_plot, x ='value', y='index', hue = 'variable', orient = 'h', ax = ax[2], hue_order = ['1.0','0.001','1e-05','1e-08'], palette=sad_colors)
		
		ax[0].set_yticklabels([self.label_dict[trait] for trait in self.summarydf_direct.index.tolist()])
		ax[0].set_ylabel('')
		ax[0].set_xlabel('Number of significant\nPC-wise components in\nthe first 20 PCs')
		ax[0].legend([],[],frameon=False)
		ax[0].set_title('Direct variance')

		ax[1].set_yticklabels(['' for trait in self.summarydf_direct.index.tolist()])
		ax[1].set_ylabel('')
		ax[1].set_xlabel('Number of significant\nPC-wise components in\nthe first 20 PCs')
		ax[1].legend([],[],frameon=False)
		ax[1].set_title('Direct-SAD covariance')

		ax[2].set_yticklabels(['' for trait in self.summarydf_direct.index.tolist()])
		ax[2].set_ylabel('')
		ax[2].set_xlabel('Number of significant\nPC-wise components in\nthe first 20 PCs')
		one_patch = mpatches.Patch(color=legend_palette['1.0'], label=r'$1$')
		three_patch = mpatches.Patch(color=legend_palette['0.001'], label=r'$0.001$')
		five_patch = mpatches.Patch(color=legend_palette['1e-05'], label=r'$1\times 10^{-5}$')
		eight_patch = mpatches.Patch(color=legend_palette['1e-08'], label=r'$1\times 10^{-8}$')
		ax[2].legend(handles=[eight_patch, five_patch, three_patch, one_patch], title = 'Ascertainment p-value')
		ax[2].set_title('SAD variance')

		for i in range(1, len(self.sps_traits)+1, 2):
  			ax[0].axhspan(i-0.5, i-1.5, facecolor='grey', alpha=0.5,zorder = 0)
  			ax[1].axhspan(i-0.5, i-1.5, facecolor='grey', alpha=0.5,zorder = 0)
  			ax[2].axhspan(i-0.5, i-1.5, facecolor='grey', alpha=0.5,zorder = 0)
  		
		ax[0].set_ylim([-0.5,16.5])
		ax[1].set_ylim([-0.5,16.5])
		ax[2].set_ylim([-0.5,16.5])
		plt.suptitle('PC-wise results for analysis: ' + label + ' ' + analysis)
		sns.despine()
		plt.tight_layout()
		fig.subplots_adjust(top=0.9)
		plt.savefig('../figures/component_significance_tables/hist.' + label + '.' + analysis + '.pdf')
		plt.clf()

	def plot_wc(self):
		counter = 0
		for analysis in [ '1kg.eur','1kg.all',]:
			for thresh in ['1.0','0.001','1e-05','1e-08']:
				outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
				for trait in self.sps_traits:
					df = pd.read_csv('../cache/component_inputs/wc/plink.wc.' + analysis + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
					sads = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
					# print(df.loc['sad_vc_pvals'][:10])
					# if analysis == '1kg.eur' and trait == 'birth_weight':
					# 	sys.exit()
					outdf.loc[trait] = outdf.loc[trait] + ',' + sads
					covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + covar
					self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
				outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
				outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
				outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
				outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
				outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
				outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))
				if counter == 0:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.xlsx') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)
					counter += 1
				else:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.xlsx', engine="openpyxl", mode = 'a') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)

			self.plot_summaries(analysis,'wc')

	def plot_nopcs(self):
		counter = 0
		for analysis in ['1kg.all', '1kg.eur']:
			for thresh in ['1.0','0.001','1e-05','1e-08']:
				outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
				for trait in self.sps_traits:
					df = pd.read_csv('../cache/component_inputs/nopcs/plink.wc.nopcs.' + analysis + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
					sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + sads
					covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + covar
					self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
				outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
				outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
				outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
				outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
				outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
				outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))

				if counter == 0:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.nopcs.xlsx') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)
					counter += 1
				else:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.nopcs.xlsx', engine="openpyxl", mode = 'a') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)

			
			self.plot_summaries(analysis,'nopcs')


	def plot_bolt(self):
		for grm in ['wpcs','nopcs']:
			counter = 0
			for analysis in ['1kg.all', '1kg.eur']:
				for thresh in ['1.0','0.001','1e-05','1e-08']:
				
					outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
					for trait in self.sps_traits:
						df = pd.read_csv('../cache/component_inputs/bolt/bolt.' + grm + '.' + analysis + '.' + trait + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
						outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
						sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
						outdf.loc[trait] = outdf.loc[trait] + ',' + sads
						covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
						outdf.loc[trait] = outdf.loc[trait] + ',' + covar
						self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
						self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
						self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
					outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
					outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
					outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
					outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
					outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
					outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))

					if counter == 0:
						with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.bolt.' + grm + '.xlsx') as writer:
							outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)
						counter += 1
					else:
						with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.bolt.' + grm + '.xlsx', engine="openpyxl", mode = 'a') as writer:
							outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)

				self.plot_summaries(analysis,'bolt.'+grm)

	def plot_ukb_and_1kg(self):
		for cohort,label in zip(['1kg.all', '1kg.eur'],['ukb.and.1kg.all.pcs','ukb.and.1kg.eur.pcs']):
			counter = 0
			for thresh in ['1.0','0.001','1e-05','1e-08']:
				outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
				for trait in self.sps_traits:
					df = pd.read_csv('../cache/component_inputs/ukb.and.1kg.pcs/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
					sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + sads
					covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + covar
					self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
				outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
				outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
				outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
				outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
				outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
				outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))
				if counter == 0:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.' + cohort + '.' + label + '.xlsx') as writer:
						outdf.to_excel(writer, sheet_name = cohort + '.' + thresh)
						counter += 1
				else:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.' + cohort + '.' + label + '.xlsx', engine="openpyxl", mode = 'a') as writer:
						outdf.to_excel(writer, sheet_name = cohort + '.' + thresh)
			self.plot_summaries(cohort,'ukb.and.1kg.pcs')

	def plot_1kg_pcs_only(self):
		for cohort,label in zip(['1kg.all', '1kg.eur'],['1kg.all.pcs.only','1kg.eur.pcs.only']):
			counter = 0
			for thresh in ['1.0','0.001','1e-05','1e-08']:
				outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
				for trait in self.sps_traits:
					df = pd.read_csv('../cache/component_inputs/1kg.pcs.only/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
					sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + sads
					covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + covar
					self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
				outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
				outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
				outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
				outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
				outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
				outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))
				if counter == 0:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.' + cohort + '.' + label + '.xlsx') as writer:
						outdf.to_excel(writer, sheet_name = cohort + '.' + thresh)
						counter += 1
				else:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.' + cohort + '.' + label + '.xlsx', engine="openpyxl", mode = 'a') as writer:
						outdf.to_excel(writer, sheet_name = cohort + '.' + thresh)
			self.plot_summaries(cohort,'pcs.only')

	def plot_maf01(self):
		counter = 0
		for analysis in ['1kg.all', '1kg.eur']:
			for thresh in ['1.0','0.001','1e-05','1e-08']:
				outdf = pd.DataFrame(np.zeros((len(self.sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = self.sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
				for trait in self.sps_traits:
					df = pd.read_csv('../cache/component_inputs/maf01/plink.wc.' + analysis + '.sps23.' + trait + '.maf01.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
					outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
					sads = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + sads
					covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
					outdf.loc[trait] = outdf.loc[trait] + ',' + covar
					self.summarydf_sad.loc[trait][thresh] = np.sum(np.where(df.loc['sad_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_direct.loc[trait][thresh] = np.sum(np.where(df.loc['direct_vc_pvals'] < 0.05,1,0)[:10])
					self.summarydf_covar.loc[trait][thresh] = np.sum(np.where(df.loc['covar_vc_pvals'] < 0.05,1,0)[:10])
				outdf['sig_sad'] = np.where(self.summarydf_sad[thresh] > 0., True, False)
				outdf['sig_direct'] = np.where(self.summarydf_direct[thresh] > 0., True, False)
				outdf['sig_covar'] = np.where(self.summarydf_covar[thresh] > 0., True, False)
				outdf.loc['total','sig_sad'] = np.sum(np.where(self.summarydf_sad[thresh] > 0., 1, 0))
				outdf.loc['total','sig_direct'] = np.sum(np.where(self.summarydf_direct[thresh] > 0., 1, 0))
				outdf.loc['total','sig_covar'] = np.sum(np.where(self.summarydf_covar[thresh] > 0., 1, 0))
				if counter == 0:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.maf01.xlsx') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)
					counter += 1
				else:
					with pd.ExcelWriter('../figures/spreadsheets/significant.pc.components.plink.wc.maf01.xlsx', engine="openpyxl", mode = 'a') as writer:
						outdf.to_excel(writer, sheet_name = analysis + '.' + thresh)

			self.plot_summaries(analysis,'wc')



