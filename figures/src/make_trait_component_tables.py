import pandas as pd 
import numpy as np 
import sys

import matplotlib.pyplot as plt


label_dict = {'pack_years_of_smoking':'Pack years of smoking',\
				'neuroticism_score':'Neuroticism score',\
				'alcohol_intake_frequency':'Alcohol intake frequency',\
				'overall_health_rating':'Overall health rating',\
				'years_of_schooling':'Years of schooling',\
				'fluid_intelligence':'Fluid intelligence',\
				'forced_vital_capacity':'Forced vital capacity',\
				'birth_weight':'Birth weight',\
				'diastolic_blood_pressure':'Diastolic blood pressure',\
				'BMI':'Body mass index',\
				'pulse_rate':'Pulse rate',\
				'waist_circumference':'Waist circumference',\
				'hand_grip_strength':'Hand grip strength',\
				'hair_color':'Hair color',\
				'hip_circumference':'Hip circumference',\
				'skin_color':'Skin color',\
				'household_income':'Household income',\
				'basal_metabolic_rate':'Basal metabolic rate',\
				'height':'Height'}
sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
		'hand_grip_strength','height','hip_circ','household_income','neuroticism_score',
		'overall_health','pack_years_smoking','pulse_rate','skin_color','waist_circ','years_schooling']

# for analysis in ['1kg.all', '1kg.eur']:
# 	for thresh in ['1.0','0.001','1e-05','1e-08']:
# 		outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
# 		for trait in sps_traits:
# 			df = pd.read_csv('../cache/component_inputs/wc/plink.wc.' + analysis + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
# 			outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 			sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 			covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 		outdf.to_csv('../figures/component_significance_tables/wc.' + analysis + '.ascp.' + thresh + '.txt',sep = '\t')
# 		outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
		
# 		for trait in sps_traits:
# 			df = pd.read_csv('../cache/component_inputs/nopcs/plink.wc.nopcs.' + analysis + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
# 			outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 			sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 			covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 		outdf.to_csv('../figures/component_significance_tables/wc.nopcs.' + analysis + '.ascp.' + thresh + '.txt',sep = '\t')

# 		for trait in sps_traits:
# 			df = pd.read_csv('../cache/component_inputs/nosingletons/plink.wc.' + analysis + '.' + trait + '.nosingletons.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
# 			outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 			sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 			covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 		outdf.to_csv('../figures/component_significance_tables/wc.nosingletons.' + analysis + '.ascp.' + thresh + '.txt',sep = '\t')

# for analysis in ['1kg.all', '1kg.eur']:
# 	for thresh in ['1.0','0.001','1e-05','1e-08']:
# 		for grm in ['wpcs','nopcs']:
# 			outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
# 			for trait in sps_traits:
# 				df = pd.read_csv('../cache/component_inputs/bolt/bolt.' + grm + '.' + analysis + '.' + trait + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
# 				outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 				sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 				outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 				covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 				outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 			outdf.to_csv('../figures/component_significance_tables/bolt.' + grm + '.' + analysis + '.ascp.' + thresh + '.txt',sep = '\t')

# for analysis in ['1kg.all', '1kg.eur']:
# 	for thresh in ['1.0','0.001','1e-05','1e-08']:
# 		for half in ['1','2']:
# 			outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
# 			for trait in sps_traits:
# 				df = pd.read_csv('../cache/component_inputs/ascertain_validate/plink.half.' + half + '.' + analysis + '.sps23.' + trait + '.aperm.1K.to.1M.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')
# 				outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 				sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 				outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 				covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 				outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 			outdf.to_csv('../figures/component_significance_tables/plink.half.' + half + '.' + analysis + '.ascp.' + thresh + '.txt',sep = '\t')

# for cohort,label in zip(['1kg.all', '1kg.eur'],['ukb.and.1kg.all.pcs','ukb.and.1kg.eur.pcs']):
# 	for thresh in ['1.0','0.001','1e-05','1e-08']:
# 		outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
# 		for trait in sps_traits:
# 			df = pd.read_csv('../cache/component_inputs/ukb.and.1kg.pcs/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')

# 			outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 			sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 			covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 		outdf.to_csv('../figures/component_significance_tables/plink.wc.' + cohort + '.' + label + '.ascp.' + thresh + '.txt',sep = '\t')

# for cohort,label in zip(['1kg.all', '1kg.eur'],['1kg.all.pcs.only','1kg.eur.pcs.only']):
# 	for thresh in ['1.0','0.001','1e-05','1e-08']:
# 		outdf = pd.DataFrame(np.zeros((len(sps_traits), len(['PC' + str(k+1) for k in range(10)]))), index = sps_traits, columns = ['PC' + str(k+1) for k in range(10)])
# 		for trait in sps_traits:
# 			df = pd.read_csv('../cache/component_inputs/1kg.pcs.only/plink.wc.' + cohort + '.sps23.' + trait + '.' + label + '.block.permutation.stats.pval.' + str(thresh) + '.txt', sep = '\t').set_index('Unnamed: 0')

# 			outdf.loc[trait] = np.where(df.loc['direct_vc_pvals'] < 0.05, 'd','')[:10]
# 			sads  = np.where(df.loc['sad_vc_pvals'] < 0.05, 's','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + sads
# 			covar = np.where(df.loc['covar_vc_pvals'] < 0.025, 'c','')[:10]
# 			outdf.loc[trait] = outdf.loc[trait] + ',' + covar
# 		outdf.to_csv('../figures/component_significance_tables/plink.wc.' + cohort + '.' + label + '.ascp.' + thresh + '.txt',sep = '\t')







