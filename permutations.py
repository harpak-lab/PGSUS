import pandas as pd 
import numpy as np 
import sys
import os
import warnings
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

class block_permutation():

	def __init__(self, block_bounds, chr_pos, chrom, standard_beta_threshed, standard_se_threshed, sib_beta_threshed, 
		sib_se_threshed, ascp, thresh, outlabel, eigenvecs, eigenvalues, emp_direct_vc, emp_sad_vc, 
		emp_covar_vc, emp_nondirect_vc, emp_standard_decomp, emp_sib_decomp, emp_diff_decomp, variance_direct_vc, 
		variance_sad_vc, variance_covar_vc, variance_nondirect_vc, outdir, pos_label, pcs_to_test=15, nperm = 1000):

		self.block_bounds= block_bounds
		self.chr_pos = chr_pos
		self.chrom = chrom
		self.standard_beta_threshed = standard_beta_threshed.reset_index(drop=True)
		self.standard_se_threshed = standard_se_threshed.reset_index(drop=True)
		self.sib_beta_threshed = sib_beta_threshed.reset_index(drop=True)
		self.sib_se_threshed = sib_se_threshed.reset_index(drop=True)
		self.nperm = nperm
		self.thresh = thresh
		self.pcs_to_test = pcs_to_test
		self.outdir = outdir
		self.outlabel = outlabel
		self.pos_label = pos_label
		self.emp_direct_vc, self.emp_sad_vc, self.emp_covar_vc, self.emp_nondirect_vc = emp_direct_vc, emp_sad_vc, emp_covar_vc, emp_nondirect_vc
		self.emp_standard_decomp, self.emp_sib_decomp, self.emp_diff_decomp = emp_standard_decomp, emp_sib_decomp, emp_diff_decomp
		self.variance_direct_vc, self.variance_sad_vc, self.variance_covar_vc, self.variance_nondirect_vc = variance_direct_vc, variance_sad_vc, variance_covar_vc, variance_nondirect_vc
		self.block_bound_process()
		self.block_dots(self.standard_beta_threshed, self.standard_se_threshed, self.sib_beta_threshed, self.sib_se_threshed, thresh, eigenvecs, self.chr_pos)
		self.permute_blocks(eigenvalues)
		self.calculate_pvalues()
		self.proportion_table()

	def block_bound_process(self):

		df = pd.read_csv(self.block_bounds, delim_whitespace = True)
		df['chr'] = df['chr'].str.replace('chr','')
		blocks = [] 
		for index, row in df.iterrows():
			snps = self.chr_pos[(self.chr_pos[self.chrom].astype(float) == float(row['chr'])) & (float(row['start']) <= self.chr_pos[self.pos_label].astype(float)) & (float(row['stop']) > self.chr_pos[self.pos_label].astype(float))].index.tolist()
			self.chr_pos.loc[snps,'block'] = int(index)

	def element_multiplier(self,x,y):
		return np.multiply(x,y)

	def summations(self,eff_sizes,eigenvecs,nsnps):
		temp = np.matmul(eff_sizes.T,eigenvecs)
		return temp		

	## First, save the relevant dot products on a per-block basis.
	## Each of these will be a matrix with #rows equal to the number of blocks and #columns
	## equal to the number of PCs. Once we have matrices like this for standard, sib, difference, and corresponding
	## errors, we can get permutation distributions by summing columns but flipping signs randomly
	## and bootstrap by summing columns that have rows randomly chosen w/ replacement from the originals.

	## Take in the same thing as the est.components, but return a list of three matrices corresponding
	## to relevant dot products per block.
	def block_dots(self, standard_beta, standard_se, sib_beta, sib_se, thresh, eigenvecs, pickrell_blocks):
		blocks = self.chr_pos['block'].unique()
		standard_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)
		sib_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)
		diff_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)
		sum_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)
		gw_err_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)
		sib_err_dots = pd.DataFrame(np.zeros((len(blocks),eigenvecs.shape[1])), index = blocks)

		for block in blocks:
			block_indices = self.chr_pos[self.chr_pos['block'] == block].index.tolist()
			standard_dots.loc[block] = self.summations(standard_beta.loc[block_indices],eigenvecs[block_indices,:], len(block_indices))
			sib_dots.loc[block] = self.summations(sib_beta.loc[block_indices],eigenvecs[block_indices,:], len(block_indices))
			diff_dots.loc[block] = self.summations(standard_beta.loc[block_indices]-sib_beta.loc[block_indices],eigenvecs[block_indices,:], len(block_indices))
			sum_dots.loc[block] = self.summations(standard_beta.loc[block_indices]+sib_beta.loc[block_indices],eigenvecs[block_indices,:], len(block_indices))
			#I think something is wrong in this decomposition
			gw_err_dots.loc[block] = self.summations(standard_se.loc[block_indices]**2,eigenvecs[block_indices,:]**2, len(block_indices))
			sib_err_dots.loc[block] = self.summations(sib_se.loc[block_indices]**2,eigenvecs[block_indices,:]**2, len(block_indices))

		self.standard_dots, self.sib_dots, self.diff_dots, self.sum_dots, self.gw_err_dots, self.sib_err_dots = standard_dots, sib_dots, diff_dots, sum_dots, gw_err_dots, sib_err_dots

	def var_comps_from_block_matrices_perm(self, standard_dots, sib_dots, diff_dots, gw_err_dots, sib_err_dots, eigenvalues):

		standard_decomp = eigenvalues * standard_dots.sum(axis=0)**2
		sib_decomp = eigenvalues * sib_dots.sum(axis=0)**2
		diff_decomp = eigenvalues * diff_dots.sum(axis=0)**2

		standard_projection = eigenvalues * standard_dots.sum(axis=0)
		sib_projection = eigenvalues * sib_dots.sum(axis=0)
		diff_projection = eigenvalues * diff_dots.sum(axis=0)

		standard_se_decomp = eigenvalues * gw_err_dots.sum(axis=0)
		sib_se_decomp = eigenvalues * sib_err_dots.sum(axis=0)

		direct_vc = sib_decomp - sib_se_decomp
		sad_vc = diff_decomp - sib_se_decomp - standard_se_decomp
		covar_vc = standard_decomp - diff_decomp - sib_decomp + 2*sib_se_decomp
		nondirect_vc = (standard_decomp - standard_se_decomp - (sib_decomp-sib_se_decomp))

		return [direct_vc, sad_vc, covar_vc, standard_decomp, sib_decomp, diff_decomp, standard_se_decomp, sib_se_decomp, standard_projection, sib_projection, diff_projection, nondirect_vc]

	def permute_blocks(self, eigenvalues):

		direct_vc_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		sad_vc_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		covar_vc_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		nondirect_vc_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))

		standard_proj_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		sib_proj_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		diff_proj_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))

		standard_decomp_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		sib_decomp_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))
		diff_decomp_perm = pd.DataFrame(np.zeros((self.nperm,len(eigenvalues))))

		for perm in range(self.nperm):

			random_signs = 2 * np.random.binomial(1,0.5,size = len(self.chr_pos['block'].unique().tolist())) - 1
			mat_standard_perm = np.apply_along_axis(self.element_multiplier, 0 ,self.standard_dots, random_signs)
			mat_sib_perm = np.apply_along_axis(self.element_multiplier, 0 ,self.sib_dots, random_signs)
			mat_diff_perm = np.apply_along_axis(self.element_multiplier, 0 ,self.diff_dots, random_signs)

			vcs_perm = self.var_comps_from_block_matrices_perm(mat_standard_perm, mat_sib_perm, mat_diff_perm, self.gw_err_dots, self.sib_err_dots, eigenvalues)
			
			direct_vc_perm.loc[perm] = vcs_perm[0]
			sad_vc_perm.loc[perm] = vcs_perm[1]
			covar_vc_perm.loc[perm] = vcs_perm[2]
			nondirect_vc_perm.loc[perm] = vcs_perm[11]

			standard_decomp_perm.loc[perm] = vcs_perm[3]
			sib_decomp_perm.loc[perm] = vcs_perm[4]
			diff_decomp_perm.loc[perm] = vcs_perm[5]

			standard_proj_perm.loc[perm] = vcs_perm[8]
			sib_proj_perm.loc[perm] = vcs_perm[9]
			diff_proj_perm.loc[perm] = vcs_perm[10]

		self.direct_vc_perm = direct_vc_perm
		self.sad_vc_perm = sad_vc_perm
		self.covar_vc_perm = covar_vc_perm
		self.nondirect_vc_perm = nondirect_vc_perm

		self.standard_proj_perm = standard_proj_perm
		self.sib_proj_perm = sib_proj_perm
		self.diff_proj_perm = diff_proj_perm

		self.standard_decomp_perm = standard_decomp_perm
		self.sib_decomp_perm = sib_decomp_perm
		self.diff_decomp_perm = diff_decomp_perm

	def calculate_pvalues(self):

		pvals_direct = np.ones(self.pcs_to_test)
		upper95_perm_direct = np.ones(self.pcs_to_test)
		lower0_perm_direct = np.ones(self.pcs_to_test)

		pvals_sad = np.ones(self.pcs_to_test)
		upper95_perm_sad = np.ones(self.pcs_to_test)
		lower0_perm_sad = np.ones(self.pcs_to_test)

		pvals_covar = np.ones(self.pcs_to_test)
		upper975_perm_covar = np.ones(self.pcs_to_test)
		lower025_perm_covar = np.ones(self.pcs_to_test)

		pvals_nondirect = np.ones(self.pcs_to_test)
		upper975_perm_nondirect = np.ones(self.pcs_to_test)
		lower025_perm_nondirect = np.ones(self.pcs_to_test)

		for k in range(self.pcs_to_test):
			
			ranked_direct = (self.direct_vc_perm[k]/np.sum(self.standard_decomp_perm, axis = 1)).sort_values().reset_index(drop=True)
			upper95_perm_direct[k] = ranked_direct.loc[949]
			lower0_perm_direct[k] = ranked_direct.loc[0]

			ranked_sad = (self.sad_vc_perm[k]/np.sum(self.standard_decomp_perm, axis = 1)).sort_values().reset_index(drop=True)
			upper95_perm_sad[k] = ranked_sad.loc[949]
			lower0_perm_sad[k] = ranked_sad.loc[0]

			ranked_covar = (self.covar_vc_perm[k]/np.sum(self.standard_decomp_perm, axis = 1)).sort_values().reset_index(drop=True)
			upper975_perm_covar[k] = ranked_covar.loc[975]
			lower025_perm_covar[k] = ranked_covar.loc[25]

			ranked_nondirect = (self.nondirect_vc_perm[k]/np.sum(self.standard_decomp_perm, axis = 1)).sort_values().reset_index(drop=True)
			upper975_perm_nondirect[k] = ranked_nondirect.loc[975]
			lower025_perm_nondirect[k] = ranked_nondirect.loc[25]

			pvals_direct[k] = np.mean(np.where((self.emp_direct_vc[k]/np.sum(self.emp_standard_decomp)) <= np.array(ranked_direct), 1, 0))
			pvals_sad[k] = np.mean(np.where((self.emp_sad_vc[k]/np.sum(self.emp_standard_decomp)) <= np.array(ranked_sad), 1, 0))
			pvals_covar[k] = np.min([float(np.mean(np.where((self.emp_covar_vc[k]/np.sum(self.emp_standard_decomp)) <= np.array(ranked_covar), 1, 0))),float(np.mean(np.where((self.emp_covar_vc[k]/np.sum(self.emp_standard_decomp)) >= np.array(ranked_covar), 1, 0)))])
			pvals_nondirect[k] = np.min([float(np.mean(np.where((self.emp_nondirect_vc[k]/np.sum(self.emp_standard_decomp)) <= np.array(ranked_nondirect), 1, 0))),float(np.mean(np.where((self.emp_nondirect_vc[k]/np.sum(self.emp_standard_decomp)) >= np.array(ranked_nondirect), 1, 0)))])

		self.pvals_direct, self.pvals_sad, self.pvals_covar, self.pvals_nondirect, self.upper95_perm_direct, self.upper95_perm_sad, self.upper975_perm_covar, self.upper975_perm_nondirect, self.lower0_perm_direct, self.lower0_perm_sad, self.lower025_perm_covar, self.lower025_perm_nondirect = pvals_direct, pvals_sad, pvals_covar, pvals_nondirect, upper95_perm_direct, upper95_perm_sad, upper975_perm_covar, upper975_perm_nondirect, lower0_perm_direct, lower0_perm_sad, lower025_perm_covar, lower025_perm_nondirect
		
	def proportion_table(self):
		
		direct_prop_by_pcs = self.emp_direct_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		direct_prop_se = self.variance_direct_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		direct_prop_se = np.sqrt(np.array(direct_prop_se).astype(float))
		direct_pvals = self.pvals_direct[:self.pcs_to_test]

		sad_prop_by_pcs = self.emp_sad_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		sad_prop_se = self.variance_sad_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		sad_prop_se = np.sqrt(np.array(sad_prop_se).astype(float))
		sad_pvals = self.pvals_sad[:self.pcs_to_test]

		covar_prop_by_pcs = self.emp_covar_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		covar_prop_se = self.variance_covar_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		covar_prop_se = np.sqrt(np.array(covar_prop_se).astype(float))
		covar_pvals = self.pvals_covar[:self.pcs_to_test]

		nondirect_prop_by_pcs = self.emp_nondirect_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		nondirect_prop_se = self.variance_nondirect_vc[:self.pcs_to_test]/np.sum(self.emp_standard_decomp)
		nondirect_prop_se = np.sqrt(np.array(nondirect_prop_se).astype(float))
		nondirect_pvals = self.pvals_nondirect[:self.pcs_to_test]

		outarray = np.vstack((direct_prop_by_pcs, direct_prop_se, direct_pvals, self.upper95_perm_direct,self.lower0_perm_direct, sad_prop_by_pcs, sad_prop_se, sad_pvals, self.upper95_perm_sad,self.lower0_perm_sad, covar_prop_by_pcs, covar_prop_se, covar_pvals, self.upper975_perm_covar,self.lower025_perm_covar, nondirect_prop_by_pcs, nondirect_prop_se, nondirect_pvals, self.upper975_perm_nondirect, self.lower025_perm_nondirect))
		outarray = np.round(outarray.astype(float),4)
		indices = ['direct_vc_estimate', 'direct_vc_estimate_se', 'direct_vc_pvals','upper95_perm_direct','lower0_perm_direct','sad_vc_estimate', 'sad_vc_estimate_se', 'sad_vc_pvals','upper95_perm_sad','lower0_perm_sad','covar_vc_estimate', 'covar_vc_estimate_se', 'covar_vc_pvals','upper975_perm_covar','lower025_perm_covar','nondirect_vc_estimate', 'nondirect_vc_estimate_se', 'nondirect_vc_pvals','upper975_perm_nondirect','lower025_perm_nondirect']
		tests_out = pd.DataFrame(outarray, index = indices, columns = ['PC' + str(i+1) for i in range(outarray.shape[1])])

		if self.outlabel == '':
			tests_out.to_csv(self.outdir + '/block.permutation.stats.pval.' + str(self.thresh) + '.txt', sep = '\t')
		else:
			tests_out.to_csv(self.outdir + '/' + self.outlabel + '.block.permutation.stats.pval.' + str(self.thresh) + '.txt', sep = '\t')
