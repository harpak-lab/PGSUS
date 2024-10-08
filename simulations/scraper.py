import pandas as pd
import numpy as np 
import argparse
import sys
import glob
import pickle 
import matplotlib.pyplot as plt 
import seaborn as sns
parser = argparse.ArgumentParser()
parser.add_argument('-a','--analysis',type=str,help='analysis to scrape', dest = 'analysis')
args = parser.parse_args()
analyses = [x for x in args.analysis.split(",")]

if 'null' in analyses:
	for architecture in ['maf','neutral']:
		env_covs = ['0.25','0.5','0.75']
		for env_var in ['0.1','0.2','0.5','0.8']:
			alpha_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
			alpha_se_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
			for env_cov in env_covs:
					#sim.neutral.architecture.pval.1.0.npcs.20.env.var.0.8.env.cov.0.75.nsib.10000.nstandard.10000.ntraitloci.200
					files = glob.glob('simulation_cache/sim.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.env.cov.' + env_cov + '.nsib.10000.nstandard.25000.ntraitloci.200/*sad.pkl')
					# print(len(files),architecture,env_cov,env_var)
					for e,file in enumerate(files):
						with open(file, 'rb') as f:
							data = pickle.load(f)
						alpha_df.loc[e,env_cov] = data['alpha']
						alpha_se_df.loc[e,env_cov] = data['alpha_se']
			alpha_df = alpha_df.iloc[:100]
			alpha_se_df = alpha_se_df.iloc[:100]
			alpha_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.alpha.txt',sep = '\t')
			alpha_se_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.alpha.se.txt',sep = '\t')

if 'thresholds' in analyses:
	for architecture in ['maf','neutral']:
		thresholds = ['1.0','0.001','0.00001','0.00000001']
		env_vars = ['0.1','0.2','0.5','0.8']
		for threshold in thresholds:
			alpha_df = pd.DataFrame(np.zeros((100,len(env_vars))),columns = env_vars)
			alpha_se_df = pd.DataFrame(np.zeros((100,len(env_vars))),columns = env_vars)
			for env_var in env_vars:
				files = glob.glob('simulation_cache/sim.' + architecture + '.architecture.pval.' + threshold + '.npcs.20.env.var.' + env_var + '.env.cov.0.5.nsib.10000.nstandard.25000.ntraitloci.200/*sad.pkl')
				for e,file in enumerate(files):
					with open(file, 'rb') as f:
						data = pickle.load(f)
					alpha_df.loc[e,env_var] = data['alpha']
					alpha_se_df.loc[e,env_var] = data['alpha_se']
			alpha_df = alpha_df.iloc[:100]
			alpha_se_df = alpha_se_df.iloc[:100]


			alpha_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.' + threshold + '.alpha.txt',sep = '\t')
			alpha_se_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.env.var.threshold.' + threshold + '.alpha.se.txt',sep = '\t')


if 'allloci' in analyses:
	for architecture in ['maf','neutral']:
		env_vars = ['0.01','0.1','0.2','0.5','0.8']
		alpha_df = pd.DataFrame(np.zeros((100,len(env_vars))),columns = env_vars)
		alpha_se_df = pd.DataFrame(np.zeros((100,len(env_vars))),columns = env_vars)
		for env_var in env_vars:
			files = glob.glob('simulation_cache/sim.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400/*.sad.pkl')
			for e,file in enumerate(files):
				with open(file, 'rb') as f:
					data = pickle.load(f)
				alpha_df.loc[e,env_var] = data['alpha']
				alpha_se_df.loc[e,env_var] = data['alpha_se']

		alpha_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.env.var.' + env_var + '.nsib.20000.nstandard.50000.ntraitloci.400.alpha.txt',sep = '\t')
		alpha_se_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.env.var.' + env_var + '.nsib.20000.nstandard.50000.ntraitloci.400.alpha.se.txt',sep = '\t')
		print(alpha_df)

if 'compare' in analyses:
	env_vars = ['0.01','0.1','0.2','0.5','0.8']
	for architecture in ['neutral','maf']:
		correlation_df = pd.DataFrame(np.zeros((1,len(env_vars))),columns = [env_vars])
		alpha_df = pd.DataFrame(np.zeros((1,len(env_vars))),columns = [env_vars])

		for env_var in env_vars:
			sibling_error = []
			standard_error = []
			files = glob.glob('simulation_cache/sim.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400/*.sad.pkl')
			for e,file in enumerate(files):
				with open(file, 'rb') as f:
					data = pickle.load(f)
				correlation_df.loc[e,env_var] = np.corrcoef((data['decomp_gwas']-data['decomp_gwas_se'])[100:],(data['decomp_sib']-data['decomp_sib_se'])[100:])[0,1]
				alpha_df.loc[e,env_var] = data['alpha']

				sib_estimates = (data['decomp_sib']-data['decomp_sib_se'])[100:]
				sib_estimands = (data['decomp_true_sib']-data['decomp_true_sib_se'])[100:]
				sibling_error.append(sib_estimands - sib_estimates)

				standard_estimates = (data['decomp_gwas']-data['decomp_gwas_se'])[100:]
				standard_estimands = (data['decomp_true_gwas']-data['decomp_true_gwas_se'])[100:]
				standard_error.append(standard_estimands - standard_estimates)
			
			standard_error = np.array(standard_error).mean(axis=0)
			sibling_error = np.array(sibling_error).mean(axis=0)
			df = pd.DataFrame(np.vstack((standard_error,sibling_error)), index = ['standard','sib'], columns = [j+101 for j in range(300)])
			df.to_csv('../figures/cache/simulation_output_cache/decomp.error.architecture.' + architecture + '.env.var.' + env_var + '.400.loci.txt',sep = '\t')
		
		mean_alpha = alpha_df.mean(axis=0)
		mean_alpha.to_csv('../figures/cache/simulation_output_cache/mean.alpha.architecture.' + architecture + '.400.loci.txt',sep = '\t') 





# if 'npcs' in analyses:
# 	for architecture in ['maf','null']:
# 		env_covs = ['0.25','0.5','0.75','1']
# 		for env_var in ['0.1','0.2','0.5','0.8']:
# 			alpha_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
# 			alpha_se_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
# 			for env_cov in env_covs:
# 					files = glob.glob('simulation_cache/sib.nosad.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.env.cov.' + env_cov + '/*.sad.pkl')
# 					for e,file in enumerate(files):
# 						with open(file, 'rb') as f:
# 							data = pickle.load(f)
# 						alpha_df.loc[e,env_cov] = data['alpha']
# 						alpha_se_df.loc[e,env_cov] = data['alpha_se']

# 			alpha_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.npcs.20.env.var.' + env_var + '.alpha.txt',sep = '\t')
# 			alpha_se_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.npcs.20.env.var.' + env_var + '.alpha.se.txt',sep = '\t')

# if 'nsib' in analyses:
# 	for architecture in ['maf','null']:
# 		env_covs = ['0.25','0.5','0.75','1']
# 		for env_var in ['0.1','0.2','0.5','0.8']:
# 			alpha_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
# 			alpha_se_df = pd.DataFrame(np.zeros((100,len(env_covs))),columns = env_covs)
# 			for env_cov in env_covs:
# 					files = glob.glob('simulation_cache/sib.nosad.' + architecture + '.architecture.pval.1.0.npcs.20.env.var.' + env_var + '.env.cov.' + env_cov + '.nsib.10000.nstandard.20000.ntraitloci.400/*.sad.pkl')
# 					print(architecture, env_cov,env_var,len(files))
# 					for e,file in enumerate(files):
# 						with open(file, 'rb') as f:
# 							data = pickle.load(f)
# 						alpha_df.loc[e,env_cov] = data['alpha']
# 						alpha_se_df.loc[e,env_cov] = data['alpha_se']

# 			alpha_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.npcs.20.env.var.' + env_var + '.nsib.10000.nstandard.20000.ntraitloci.400.alpha.txt',sep = '\t')
# 			alpha_se_df.to_csv('simulation_output_cache/null.' + architecture + '.architecture.npcs.20.env.var.' + env_var + '.nsib.10000.nstandard.20000.ntraitloci.400.alpha.se.txt',sep = '\t')




# 	env_vars = ['0.1','0.5']
# 	standard_n = ['2000','10000','20000','50000','100000','200000']
# 	for variance in env_vars:
# 		alpha_df = pd.DataFrame(np.zeros((100,len(standard_n))),columns = standard_n)
# 		ols_alpha_df = pd.DataFrame(np.zeros((100,len(standard_n))),columns = standard_n)
# 		for n in standard_n:
# 			alpha_vc_df = pd.DataFrame(np.zeros((2,300)), index = ['sib','standard'], columns = [i for i in range(100,400)])
# 			alpha_sibling_vcs_df = pd.DataFrame(np.zeros((100,300)),columns = [i for i in range(100,400)])
# 			alpha_standard_vcs_df = pd.DataFrame(np.zeros((100,300)),columns = [i for i in range(100,400)])

# 			files = glob.glob('simulation_cache/sib.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0/*.sad.pkl')[:100]
# 			print(files)
# 			with open(files[0], 'rb') as f:
# 				data = pickle.load(f)			
			
# 			npcs = len(data['direct_variance_component_error'])
# 			nsims = len(files)

# 			total_var_df = pd.DataFrame(np.zeros((nsims,3)), columns = ['direct','sad','covar'])
# 			direct_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))
# 			sad_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))
# 			covar_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))

# 			for e,file in enumerate(files):
# 				with open(file, 'rb') as f:
# 					data = pickle.load(f)
				
# 				alpha_df.loc[e,n] = data['alpha']
# 				ols_alpha_df.loc[e,n] = data['alpha_ols']
# 				total_var_df.loc[e] = [np.sum(data['direct_vc'][100:]),np.sum(data['sad_vc'][100:]),np.sum(data['covar_vc'][100:])]
# 				alpha_sibling_vcs_df.loc[e] = data['lmdf']['sib_vc']
# 				alpha_standard_vcs_df.loc[e] = data['lmdf']['standard_vc']
# 				direct_vc_error_df.loc[e] = data['direct_variance_component_error']
# 				sad_vc_error_df.loc[e] = data['sad_variance_component_error']
# 				covar_vc_error_df.loc[e] = data['covar_variance_component_error']

# 			alpha_vc_df.loc['sib'] = alpha_sibling_vcs_df.mean(axis=0)
# 			alpha_vc_df.loc['standard'] = alpha_standard_vcs_df.mean(axis=0)
# 			alpha_vc_df = alpha_vc_df.T

# 			alpha_vc_df.to_csv('fig_cache/nratios/lmdf.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)			

# 			out_component = pd.DataFrame(np.zeros((npcs,6)), columns = ['direct','direct_se','sad','sad_se','covar','covar_se',])
# 			out_component['direct'] = direct_vc_error_df.mean(axis = 0)
# 			out_component['direct_se'] = direct_vc_error_df.sem(axis = 0)
# 			out_component['sad'] = sad_vc_error_df.mean(axis = 0)
# 			out_component['sad_se'] = sad_vc_error_df.sem(axis = 0)
# 			out_component['covar'] = covar_vc_error_df.mean(axis = 0)
# 			out_component['covar_se'] = covar_vc_error_df.sem(axis = 0)

# 			# alpha_vc_df.to_csv('fig_cache/nratios/lmdf.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)
# 			out_component.to_csv('fig_cache/nratios/component.error.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)
# 			total_var_df.to_csv('fig_cache/nratios/total.var.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)

# 		alpha_df.to_csv('fig_cache/nratios/alpha.estimates.env.var.' + variance + '.txt',sep = '\t',index = False)
# 		ols_alpha_df.to_csv('fig_cache/nratios/ols.alpha.estimates.env.var.' + variance + '.txt',sep = '\t',index = False)

# 	for variance in env_vars:
# 		mse_df = pd.DataFrame(np.zeros((100,len(standard_n))),index = [i for i in range(100)], columns = standard_n)
# 		beta_se_diff_df = pd.DataFrame(np.zeros((100,len(standard_n))),index = [i for i in range(100)], columns = standard_n)
# 		beta_true2_df = pd.DataFrame(np.zeros((100,len(standard_n))),index = [i for i in range(100)], columns = standard_n)
# 		beta2_df = pd.DataFrame(np.zeros((100,len(standard_n))),index = [i for i in range(100)], columns = standard_n)
# 		sib_beta_se_diff_df = pd.DataFrame(np.zeros((100,len(standard_n))),index = [i for i in range(100)], columns = standard_n)
		
# 		for n in standard_n:
# 			files = glob.glob('simulation_cache/sib.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0/*.gwas.pkl')[:100]

# 			for e,file in enumerate(files):
# 				with open(file, 'rb') as f:
# 					data = pickle.load(f)

# 				mse_df.loc[e,n] = np.mean((data['true_standard_beta'] - data['gwas_combined']['beta'])**2)
# 				beta_se_diff_df.loc[e,n] = np.mean((data['gwas_combined']['beta']**2)-(data['gwas_combined']['se']**2))
# 				beta_true2_df.loc[e,n] = np.mean(data['true_standard_beta']**2)
# 				beta2_df.loc[e,n] = np.mean(data['gwas_combined']['beta']**2)
# 				sib_beta_se_diff_df.loc[e,n] = np.mean((data['sibgwas_uncorr']['beta']**2)-(data['sibgwas_uncorr']['se']**2))

# 		ratiodf = beta_se_diff_df/sib_beta_se_diff_df

# 		true_diff_df = beta_true2_df - beta_se_diff_df
# 		diff_df = beta2_df - beta_se_diff_df
# 		ratiodf.to_csv('fig_cache/nratios/ratio.nstandard.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)
# 		mse_df.to_csv('fig_cache/nratios/mse.nstandard.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)
# 		true_diff_df.to_csv('fig_cache/nratios/beta.2.se.2.diff.nstandard.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)

# if 'nloci' in analyses:
# 	env_vars = ['0.1','0.5']
# 	nloci = ['1','2','4','20','40','100','200','300']
# 	for variance in env_vars:
# 		alpha_df = pd.DataFrame(np.zeros((100,len(nloci))),columns = nloci)
# 		ols_alpha_df = pd.DataFrame(np.zeros((100,len(nloci))),columns = nloci)
# 		for n in nloci:
# 			alpha_vc_df = pd.DataFrame(np.zeros((2,300)), index = ['sib','standard'], columns = [i for i in range(100,400)])
# 			alpha_sibling_vcs_df = pd.DataFrame(np.zeros((100,300)),columns = [i for i in range(100,400)])
# 			alpha_standard_vcs_df = pd.DataFrame(np.zeros((100,300)),columns = [i for i in range(100,400)])
			
# 			files = glob.glob('simulation_cache/nloci.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0/*.sad.pkl')[:100]
# 			with open(files[0], 'rb') as f:
# 				data = pickle.load(f)			
			
# 			npcs = len(data['direct_variance_component_error'])
# 			nsims = len(files)

# 			total_var_df = pd.DataFrame(np.zeros((nsims,3)), columns = ['direct','sad','covar'])
# 			direct_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))
# 			sad_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))
# 			covar_vc_error_df = pd.DataFrame(np.zeros((nsims,npcs)))

# 			for e,file in enumerate(files):
# 				with open(file, 'rb') as f:
# 					data = pickle.load(f)
				
# 				alpha_df.loc[e,n] = data['alpha']
# 				ols_alpha_df.loc[e,n] = data['alpha_ols']
# 				total_var_df.loc[e] = [np.sum(data['direct_vc'][100:]),np.sum(data['sad_vc'][100:]),np.sum(data['covar_vc'][100:])]
# 				alpha_sibling_vcs_df.loc[e] = data['lmdf']['sib_vc']
# 				alpha_standard_vcs_df.loc[e] = data['lmdf']['standard_vc']
# 				direct_vc_error_df.loc[e] = data['direct_variance_component_error']
# 				sad_vc_error_df.loc[e] = data['sad_variance_component_error']
# 				covar_vc_error_df.loc[e] = data['covar_variance_component_error']

# 			alpha_vc_df.loc['sib'] = alpha_sibling_vcs_df.mean(axis=0)
# 			alpha_vc_df.loc['standard'] = alpha_standard_vcs_df.mean(axis=0)
# 			alpha_vc_df = alpha_vc_df.T

# 			alpha_vc_df.to_csv('fig_cache/nloci/lmdf.nloci.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)			

# 			out_component = pd.DataFrame(np.zeros((npcs,6)), columns = ['direct','direct_se','sad','sad_se','covar','covar_se',])
# 			out_component['direct'] = direct_vc_error_df.mean(axis = 0)
# 			out_component['direct_se'] = direct_vc_error_df.sem(axis = 0)
# 			out_component['sad'] = sad_vc_error_df.mean(axis = 0)
# 			out_component['sad_se'] = sad_vc_error_df.sem(axis = 0)
# 			out_component['covar'] = covar_vc_error_df.mean(axis = 0)
# 			out_component['covar_se'] = covar_vc_error_df.sem(axis = 0)

# 			# alpha_vc_df.to_csv('fig_cache/nratios/lmdf.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)
# 			out_component.to_csv('fig_cache/nloci/component.error.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t', index = False)
# 			total_var_df.to_csv('fig_cache/nloci/total.var.nstandard.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)

# 		alpha_df.to_csv('fig_cache/nloci/alpha.estimates.env.var.' + variance + '.txt',sep = '\t',index = False)
# 		ols_alpha_df.to_csv('fig_cache/nloci/ols.alpha.estimates.env.var.' + variance + '.txt',sep = '\t',index = False)

# 	for variance in env_vars:
# 		mse_df = pd.DataFrame(np.zeros((100,len(nloci))),index = [i for i in range(100)], columns = nloci)
# 		beta_se_diff_df = pd.DataFrame(np.zeros((100,len(nloci))),index = [i for i in range(100)], columns = nloci)
# 		beta_true2_df = pd.DataFrame(np.zeros((100,len(nloci))),index = [i for i in range(100)], columns = nloci)
# 		beta2_df = pd.DataFrame(np.zeros((100,len(nloci))),index = [i for i in range(100)], columns = nloci)
# 		sib_beta_se_diff_df = pd.DataFrame(np.zeros((100,len(nloci))),index = [i for i in range(100)], columns = nloci)
		
# 		for n in nloci:
# 			files = glob.glob('simulation_cache/nloci.' + n + '.env.var.' + variance + '.env.covar.0.5.pval.1.0/*.gwas.pkl')[:100]

# 			for e,file in enumerate(files):
# 				with open(file, 'rb') as f:
# 					data = pickle.load(f)

# 				mse_df.loc[e,n] = np.mean((data['true_standard_beta'] - data['gwas_combined']['beta'])**2)
# 				beta_se_diff_df.loc[e,n] = np.mean((data['gwas_combined']['beta']**2)-(data['gwas_combined']['se']**2))
# 				beta_true2_df.loc[e,n] = np.mean(data['true_standard_beta']**2)
# 				beta2_df.loc[e,n] = np.mean(data['gwas_combined']['beta']**2)
# 				sib_beta_se_diff_df.loc[e,n] = np.mean((data['sibgwas_uncorr']['beta']**2)-(data['sibgwas_uncorr']['se']**2))

# 		ratiodf = beta_se_diff_df/sib_beta_se_diff_df

# 		true_diff_df = beta_true2_df - beta_se_diff_df
# 		diff_df = beta2_df - beta_se_diff_df
# 		ratiodf.to_csv('fig_cache/nloci/ratio.nloci.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)
# 		mse_df.to_csv('fig_cache/nloci/mse.nloci.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)
# 		true_diff_df.to_csv('fig_cache/nloci/beta.2.se.2.diff.nloci.env.var.' + variance + '.env.covar.0.5.pval.1.0',sep = '\t',index = False)



