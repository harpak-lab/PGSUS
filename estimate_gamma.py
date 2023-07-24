import numpy as np
import pandas as pd
import sys
import statsmodels.formula.api as smf
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

mass = importr('MASS')
base = importr('base')
stats = importr('stats')


class SAD(object):

	def __init__(self, gwas_data_dict,thresh, pc_lower_bound, weight = False):
		self.gwas_beta = gwas_data_dict['gwas_combined']['beta']

		self.gwas_se = gwas_data_dict['gwas_combined']['se']
		self.ascertainment_p = gwas_data_dict['gwas_combined']['p']

		self.sib_beta = gwas_data_dict['sibgwas_uncorr']['beta']
		self.sib_se = gwas_data_dict['sibgwas_uncorr']['se']

		self.thresh = thresh
		self.pc_lower_bound = pc_lower_bound

		self.eigenvectors = gwas_data_dict['prs_pcs_vectors']
		self.eigenvalues = gwas_data_dict['prs_pcs_values']
		
		self.weight = weight
		self.pc_upper_bound = self.eigenvalues.shape[0]
		self.nloci = self.gwas_beta.shape[0]
		self.direct_variance_component,self.sad_variance_component, self.covar_variance_component, self.decomp_gwas, self.decomp_sib, self.decomp_diff, self.gwas_avg_se2,self.sib_avg_se2, self.proj_gwas, self.proj_sib, self.proj_diff, self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, self.beta_sum, self.gamma, self.se = self.estimate_components_and_gamma(self.gwas_beta,self.gwas_se,self.sib_beta,self.sib_se,self.ascertainment_p, self.thresh,self.eigenvectors,self.eigenvalues,self.pc_upper_bound, weight=False, pc_lower_bound = self.pc_lower_bound)
	
	def element_multiplier(self,x,y):
		return np.multiply(x,y)
	
	#need to apply this on every column
	def prs_decomp(self,eff_sizes,eigenvalues,eigenvecs):
		eff_sizes = np.array(eff_sizes).reshape(self.nloci)
		temp = np.apply_along_axis(self.element_multiplier, 0 , eigenvecs, eff_sizes)
		temp = np.sum(temp, axis = 0)**2
		temp = np.multiply(np.array(eigenvalues).reshape(self.pc_upper_bound),temp)
		return temp

	def error_decomp(self,ses,eigenvalues,eigenvecs):
		ses = np.array(ses).reshape(self.nloci)
		temp = np.apply_along_axis(self.element_multiplier, 0 , eigenvecs**2, ses**2)
		temp = np.sum(temp, axis = 0)
		temp = np.multiply(np.array(eigenvalues).reshape(self.pc_upper_bound),temp)
		return temp

	def prs_decomp_unsq(self,eff_sizes,eigenvalues,eigenvecs):
		eff_sizes = np.array(eff_sizes).reshape(self.nloci)
		temp = np.apply_along_axis(self.element_multiplier, 0 , eigenvecs, eff_sizes)
		temp = np.sum(temp, axis = 0)
		temp = np.multiply(np.array(eigenvalues).reshape(self.pc_upper_bound),temp)
		return temp

	#zeros out the effect size estimates (beta) where p is greater than threshold
	def beta_p_thresh(self,betas,p,thresh):
		filt_p = p[p >= thresh].index.tolist()
		betas.loc[betas.index.isin(filt_p)] = 0
		return betas
		
	def se2_avg_p_thresh(self,ses,p,thresh):
		filt_p = p[p >= thresh].index.tolist()
		return np.nanmean(ses.loc[~ses.index.isin(filt_p)]**2)

	def estimate_components_and_gamma(self,gwas_beta,gwas_se,sib_beta,sib_se,ascertainment_p,thresh,eigenvecs,eigenvalues,
		pc_upper_bound, weight=False, pc_lower_bound = 100):
		
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
		# print(gwas_beta_threshed)
		#perform the decomposition for each component of interest
		decomp_gwas = self.prs_decomp(gwas_beta_threshed,eigenvalues,eigenvecs)
		decomp_sib = self.prs_decomp(sib_beta_threshed,eigenvalues,eigenvecs)
		decomp_diff = self.prs_decomp(gwas_beta_threshed - sib_beta_threshed, eigenvalues, eigenvecs)

		#perform the decomposition for each component of interest, unsquared
		proj_gwas = self.prs_decomp_unsq(gwas_beta_threshed,eigenvalues,eigenvecs)
		proj_sib = self.prs_decomp_unsq(sib_beta_threshed,eigenvalues,eigenvecs)
		proj_diff = self.prs_decomp_unsq(gwas_beta_threshed - sib_beta_threshed,eigenvalues,eigenvecs)

		#decompose the standard errors of the gwas and sibling standard errors
		decomp_gwas_se = self.error_decomp(gwas_se_threshed, eigenvalues, eigenvecs)
		decomp_sib_se = self.error_decomp(sib_se_threshed, eigenvalues, eigenvecs)

		#get the variance components for direct, sad, and their covariance
		direct_variance_component = decomp_sib - decomp_sib_se
		sad_variance_component = decomp_diff - decomp_sib_se - decomp_gwas_se
		covar_variance_component = decomp_gwas - decomp_diff - decomp_sib + 2*decomp_sib_se

		#get variance of each variance component - woah, meta
		variance_direct_vc = 2*(decomp_sib_se**2)
		variance_sad_vc = 2*((decomp_sib_se+decomp_gwas_se)**2)
		variance_covar_vc = 4*(decomp_sib_se+decomp_gwas_se) + 8*(decomp_sib_se**2)

		sib_vc = direct_variance_component
		gwas_vc = decomp_gwas - decomp_gwas_se

		startdf = np.vstack((sib_vc,gwas_vc,eigenvalues))

		lmdf = pd.DataFrame(data=startdf, index = ['sib_vc','gwas_vc','weights']).T
		lmdf = lmdf.astype(float)
		lmdf = lmdf.iloc[pc_lower_bound:pc_upper_bound]

		indices = [i for i in range(pc_lower_bound,pc_upper_bound)]

		if self.weight is False:
			pandas2ri.activate()
			robjects.globalenv['dataframe'] = lmdf
			M = stats.lm(formula='gwas_vc ~ sib_vc - 1',data = base.as_symbol('dataframe'))
			gamma = np.sqrt(base.summary(M).rx2('coefficients')[0,0])
			se = base.summary(M).rx2('coefficients')[0,1]

		else:
			pandas2ri.activate()
			robjects.globalenv['dataframe'] = lmdf
			M = stats.lm(formula='gwas_vc ~ sib_vc - 1',weights = lmdf.loc[indices]['weights'],data = base.as_symbol('dataframe'))
			gamma = np.sqrt(base.summary(M).rx2('coefficients')[0,0])
			se = base.summary(M).rx2('coefficients')[0,1]

		return direct_variance_component, sad_variance_component, covar_variance_component, decomp_gwas, decomp_sib, decomp_diff, gwas_avg_se2, sib_avg_se2, proj_gwas, proj_sib, proj_diff, variance_direct_vc, variance_sad_vc, variance_covar_vc, np.sum(gwas_beta_threshed[gwas_beta_threshed != 0]), gamma, se

	def outputs(self):
		return {'direct_vc':self.direct_variance_component,'sad_vc':self.sad_variance_component, 'covar_vc':self.covar_variance_component, 
				'decomp_gwas':self.decomp_gwas, 'decomp_sib':self.decomp_sib, 'decomp_diff':self.decomp_diff,
				'gwas_avg_se2':self.gwas_avg_se2,'sib_avg_se2':self.sib_avg_se2,
				'proj_gwas':self.proj_gwas, 'proj_sib':self.proj_sib, 'proj_diff':self.proj_diff,
				'var_direct_vc':self.variance_direct_vc, 'var_sad_vc':self.variance_sad_vc, 'var_covar_vc':self.variance_covar_vc,
				'beta_sum':self.beta_sum, 'gamma':self.gamma, 'se':self.se}


