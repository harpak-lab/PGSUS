import numpy as np
import pandas as pd
import sys
import statsmodels.formula.api as smf
from scipy.odr import *
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

class SAD(object):

	def __init__(self, gwas_data_dict,thresh, pc_lower_bound,sib_se_multiplier=1):
		self.gwas_beta = gwas_data_dict['gwas_combined']['beta']
		self.gwas_se = gwas_data_dict['gwas_combined']['se']
		self.ascertainment_p = gwas_data_dict['gwas_combined']['p']
		self.sib_beta = gwas_data_dict['sibgwas_uncorr']['beta']
		self.sib_se = gwas_data_dict['sibgwas_uncorr']['se']*sib_se_multiplier
		self.thresh = thresh
		self.pc_lower_bound = pc_lower_bound

		self.eigenvectors = gwas_data_dict['target_eigenvectors']
		self.eigenvalues = gwas_data_dict['target_eigenvalues']

		self.true_gwas_beta = pd.Series(gwas_data_dict['true_standard_beta'])
		self.true_gwas_se = pd.Series(gwas_data_dict['true_standard_ses'])
		self.true_sib_beta = pd.Series(gwas_data_dict['true_sib_beta'])
		self.true_sib_se = pd.Series(gwas_data_dict['true_sib_ses'])

		self.pc_upper_bound = self.eigenvalues.shape[0]
		self.nloci = self.gwas_beta.shape[0]
		self.direct_variance_component, self.sad_variance_component, self.covar_variance_component, self.decomp_gwas, self.decomp_sib, self.decomp_diff, self.gwas_avg_se2,self.sib_avg_se2, self.proj_gwas, self.proj_sib, self.proj_diff, self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, self.beta_sum, self.alpha, self.alpha_se = self.estimate_components_and_alpha(self.gwas_beta,self.gwas_se,self.sib_beta,self.sib_se,self.ascertainment_p, self.thresh,self.eigenvectors,self.eigenvalues,self.true_gwas_beta, self.true_gwas_se, self.true_sib_beta, self.true_sib_se, self.pc_upper_bound, pc_lower_bound = self.pc_lower_bound)

	def element_multiplier(self,x,y):
		return np.multiply(x,y)

	def prs_decomp(self,eff_sizes,eigenvalues,eigenvectors):
		eff_sizes = np.array(eff_sizes)
		temp = np.apply_along_axis(self.element_multiplier,0,eigenvectors,eff_sizes)
		temp = np.sum(temp, axis = 0)
		temp = np.power(temp,2)
		temp = np.multiply(eigenvalues,temp)
		return temp

	def error_decomp(self,ses,eigenvalues,eigenvectors):
		ses = np.power(np.array(ses).reshape(self.eigenvectors.shape[0]),2)
		temp = np.apply_along_axis(self.element_multiplier, 0, np.power(eigenvectors,2), ses)
		temp = np.sum(temp, axis = 0)
		temp = np.multiply(eigenvalues,temp)
		return temp
    
	def deming_error_decomp(self,errors,eigenvalues,eigenvecs):
		errors = np.array(errors).reshape(self.eigenvectors.shape[0])
		temp = np.sum(np.apply_along_axis(self.element_multiplier, 0, np.power(eigenvecs,2), errors),axis = 0)
		temp = np.power(temp,2)
		temp = 2*np.multiply(np.power(eigenvalues,2),temp)
		return temp

	def prs_decomp_unsq(self,eff_sizes,eigenvalues,eigenvectors):
		eff_sizes = np.array(eff_sizes).reshape(self.eigenvectors.shape[0])
		temp = np.apply_along_axis(self.element_multiplier, 0, eigenvectors,eff_sizes)
		temp = np.sum(temp, axis = 0)
		temp = np.multiply(eigenvalues,temp)
		return temp
	
	def f(self,B,x):
		return B[0]*x

	def beta_p_thresh(self,betas,p,thresh):
		filt_p = p[p >= thresh].index.tolist()
		betas.loc[betas.index.isin(filt_p)] = 0
		return betas
		
	def se2_avg_p_thresh(self,ses,p,thresh):
		filt_p = p[p >= thresh].index.tolist()
		return np.nanmean(ses.loc[~ses.index.isin(filt_p)]**2)

	def estimate_components_and_alpha(self,gwas_beta,gwas_se,sib_beta,sib_se,ascertainment_p,thresh,eigenvectors,eigenvalues,
		true_gwas_beta, true_gwas_se, true_sib_beta, true_sib_se, pc_upper_bound, pc_lower_bound = 100):
		
		#get all necessary statistics from the population gwas
		#get gwas effects that have p values less than threshold
		gwas_beta_threshed = self.beta_p_thresh(gwas_beta, ascertainment_p, thresh)
		#get gwas SEs that have p values less than threshold
		gwas_se_threshed = self.beta_p_thresh(gwas_se, ascertainment_p, thresh)
		#get gwas SEs that have p values less than threshold
		gwas_avg_se2 = self.se2_avg_p_thresh(gwas_se, ascertainment_p, thresh)
		
		#get all necessary statistics from the sib gwas
		#get sib effects that have p values less than threshold
		sib_beta_threshed = self.beta_p_thresh(sib_beta, ascertainment_p, thresh)
		#get sib effects that have p values less than threshold
		sib_se_threshed = self.beta_p_thresh(sib_se, ascertainment_p, thresh)
		#get sib SEs that have p values less than threshold
		sib_avg_se2 = self.se2_avg_p_thresh(sib_se, ascertainment_p, thresh)

		divided = np.array((sib_beta_threshed**2-sib_se_threshed**2))
		divided = divided[~np.isinf(divided)]
		divided = divided[~np.isnan(divided)]
		divided = divided[divided != 0]

		divided = np.array((gwas_beta_threshed**2-gwas_se_threshed**2))
		divided = divided[~np.isinf(divided)]
		divided = divided[~np.isnan(divided)]
		divided = divided[divided != 0]

		#perform the decomposition for each component of interest
		decomp_gwas = self.prs_decomp(gwas_beta_threshed,eigenvalues,eigenvectors)
		decomp_sib = self.prs_decomp(sib_beta_threshed,eigenvalues,eigenvectors)
		decomp_diff = self.prs_decomp(gwas_beta_threshed - sib_beta_threshed, eigenvalues, eigenvectors)

		#perform the decomposition for each component of interest, unsquared
		proj_gwas = self.prs_decomp_unsq(gwas_beta_threshed,eigenvalues,eigenvectors)
		proj_sib = self.prs_decomp_unsq(sib_beta_threshed,eigenvalues,eigenvectors)
		proj_diff = self.prs_decomp_unsq(gwas_beta_threshed - sib_beta_threshed,eigenvalues,eigenvectors)

		#decompose the standard errors of the gwas and sibling standard errors
		decomp_gwas_se = self.error_decomp(gwas_se_threshed, eigenvalues, eigenvectors)
		decomp_sib_se = self.error_decomp(sib_se_threshed, eigenvalues, eigenvectors)

		# #get the variance components for direct, sad, and their covariance
		direct_variance_component = decomp_sib - decomp_sib_se
		sad_variance_component = decomp_diff - decomp_sib_se - decomp_gwas_se
		covar_variance_component = decomp_gwas - decomp_diff - decomp_sib +2*decomp_sib_se
		standard_variance_component = decomp_gwas - decomp_gwas_se

		#get variance of each variance component - woah, meta
		variance_direct_vc = 2*(decomp_sib_se**2)
		variance_sad_vc = 2*((decomp_sib_se+decomp_gwas_se)**2)
		variance_covar_vc = 4*(decomp_sib_se+decomp_gwas_se) + 8*(decomp_sib_se**2)
		
		variance_standard_vc = np.power(2*((decomp_gwas_se)**2),0.5)
		variance_sib_vc = np.power(2*((decomp_sib_se)**2),0.5)

		indices = [i for i in range(pc_lower_bound,pc_upper_bound)]

		startdf = np.vstack((standard_variance_component,direct_variance_component,variance_standard_vc,variance_sib_vc))
		lmdf = pd.DataFrame(data=startdf, index = ['standard_vc','sib_vc','var_standard_vc','var_sib_vc']).T
		lmdf = lmdf.astype(float)
		lmdf = lmdf.iloc[pc_lower_bound:pc_upper_bound]
		
		linear = Model(self.f)
		mydata = RealData(x=lmdf['sib_vc'],y=lmdf['standard_vc'],sx=lmdf['var_sib_vc'], sy=lmdf['var_standard_vc'])
		myodr = ODR(mydata, linear, beta0 = [0.])
		myoutput = myodr.run()
		alpha = np.sqrt(np.abs(myoutput.beta[0]))
		alpha_se = np.sqrt(myoutput.sd_beta[0]/lmdf.shape[0])
		self.lmdf = lmdf

		model = smf.ols('standard_vc ~ 0 + sib_vc',data = lmdf).fit()
		self.ols_alpha = model.params[0]

		####The code below gets the estimands for each component using the true beta and corresponding SEs
		#get all necessary statistics from the population gwas
		#get gwas effects that have p values less than threshold
		true_gwas_beta_threshed = self.beta_p_thresh(true_gwas_beta, ascertainment_p, thresh)
		#get gwas SEs that have p values less than threshold
		true_gwas_se_threshed = self.beta_p_thresh(true_gwas_se, ascertainment_p, thresh)
		#get gwas SEs that have p values less than threshold
		true_gwas_avg_se2 = self.se2_avg_p_thresh(true_gwas_se, ascertainment_p, thresh)
		
		#get all necessary statistics from the sib gwas
		#get sib effects that have p values less than threshold
		true_sib_beta_threshed = self.beta_p_thresh(true_sib_beta, ascertainment_p, thresh)
		#get sib effects that have p values less than threshold
		true_sib_se_threshed = self.beta_p_thresh(true_sib_se, ascertainment_p, thresh)
		#get sib SEs that have p values less than threshold
		true_sib_avg_se2 = self.se2_avg_p_thresh(true_sib_se, ascertainment_p, thresh)
		# print(gwas_beta_threshed)
		#perform the decomposition for each component of interest
		decomp_true_gwas = self.prs_decomp(true_gwas_beta_threshed,eigenvalues,eigenvectors)
		decomp_true_sib = self.prs_decomp(true_sib_beta_threshed,eigenvalues,eigenvectors)
		decomp_true_diff = self.prs_decomp(true_gwas_beta_threshed - true_sib_beta_threshed, eigenvalues, eigenvectors)

		#decompose the standard errors of the gwas and sibling standard errors
		decomp_true_gwas_se = self.error_decomp(true_gwas_se_threshed, eigenvalues, eigenvectors)
		decomp_true_sib_se = self.error_decomp(true_sib_se_threshed, eigenvalues, eigenvectors)

		# #get the variance components for direct, sad, and their covariance
		true_direct_variance_component = decomp_true_sib - decomp_true_sib_se
		true_sad_variance_component = decomp_true_diff - decomp_true_sib_se - decomp_true_gwas_se
		true_covar_variance_component = decomp_true_gwas - decomp_true_diff - decomp_true_sib +2*decomp_true_sib_se
		true_standard_variance_component = decomp_true_gwas - decomp_true_gwas_se

		self.direct_variance_component_error = true_direct_variance_component - direct_variance_component
		self.sad_variance_component_error = true_sad_variance_component - sad_variance_component
		self.covar_variance_component_error = true_covar_variance_component - covar_variance_component

		return direct_variance_component, sad_variance_component, covar_variance_component, decomp_gwas, decomp_sib, decomp_diff, gwas_avg_se2, sib_avg_se2, proj_gwas, proj_sib, proj_diff, variance_direct_vc, variance_sad_vc, variance_covar_vc, np.sum(gwas_beta_threshed[gwas_beta_threshed != 0]), alpha, alpha_se

	def outputs(self):
		return {'direct_vc':self.direct_variance_component,'sad_vc':self.sad_variance_component, 'covar_vc':self.covar_variance_component, 
				'decomp_gwas':self.decomp_gwas, 'decomp_sib':self.decomp_sib, 'decomp_diff':self.decomp_diff,
				'gwas_avg_se2':self.gwas_avg_se2,'sib_avg_se2':self.sib_avg_se2,
				'proj_gwas':self.proj_gwas, 'proj_sib':self.proj_sib, 'proj_diff':self.proj_diff,
				'var_direct_vc':self.variance_direct_vc, 'var_sad_vc':self.variance_sad_vc, 'var_covar_vc':self.variance_covar_vc,
				'beta_sum':self.beta_sum, 'alpha':self.alpha, 'alpha_se':self.alpha_se, 'direct_variance_component_error':self.direct_variance_component_error,
				'sad_variance_component_error':self.sad_variance_component_error, 'covar_variance_component_error':self.covar_variance_component_error,
				'lmdf':self.lmdf,'alpha_ols':self.ols_alpha}


