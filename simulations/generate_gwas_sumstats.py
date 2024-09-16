import numpy as np
import pandas as pd
import sys
import statsmodels.formula.api as smf
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import scipy
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

class GWAS(object):

	def __init__(self, nloci, nanc, n_standard_gwas, n_sib_gwas, n_target_cohort, n_trait_loci, 
		direct_effect_variance, indirect_effect_variance, direct_indirect_covariance, drift_param, 
		env_variance, env_covariance, gwas_pcs, boot = False, 
		neutral_architecture = False, covs = None, int_naught = False):

		self.int_naught = int_naught
		self.boot = boot
		self.nloci = nloci
		self.n_trait_loci = n_trait_loci
		self.n_standard_gwas = n_standard_gwas
		self.n_sib_gwas = n_sib_gwas
		self.n_target_cohort = n_target_cohort
		self.gwas_pcs = gwas_pcs
		self.env_variance = env_variance
		self.env_covariance = env_covariance
		self.direct_effect_variance = direct_effect_variance
		self.indirect_effect_variance = indirect_effect_variance
		self.direct_indirect_covariance = direct_indirect_covariance
		self.gwas_combined, self.sibgwas_uncorr_diff, self.pc_sol_prs_vectors, self.pc_sol_prs_values, self.prs_gwas_comb, self.prs_sib,self.trait_prs,self.af = self.runner(self.nloci, self.n_standard_gwas, n_sib_gwas, self.n_target_cohort, 
			self.n_trait_loci, direct_effect_variance, indirect_effect_variance, 
			direct_indirect_covariance, drift_param, env_variance, env_covariance, 
			boot = boot, nanc = nanc, neutral_architecture = neutral_architecture)
		self.outputs()

	def extract_lm_coeffs(self,gens, phens, covs = None, int_naught = False, boot = False):
		#Run a GWAS for one locus.
		#Gens and phens are vectors of genotypes and phenotypes.
		#optionally, the user can add covariates to the regression,
		#set an intercept of 0 (for sib regression), or get SEs via bootstrap
		#Bootstrap is done w/ 500 bootstrap samples, and the SE is *not* used for
		#calculating the p value at this moment, just for the SE estimate.
		# gens = (gens - gens.mean())/(gens.std())
		#check to see if SNP is monomorphic


		if np.var(gens) == 0:
			 return np.array([0,np.nan,np.nan,1])

		else:
			#otherwise run regression, either without covariates or with them
			output_sumstats = []
			if int_naught == False:

				if covs.shape[1] == 0:
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
					try:
						model = smf.ols('phenotype ~ genotype', data=dataframe).fit()
						t = float(model.params[1])/model.bse[1]
						p = scipy.stats.norm.sf(abs(t))*2
						return [model.params[1], model.bse[1], t, p] 						
					except ValueError:
						return np.array([0,np.nan,np.nan,1])

				else:
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
					for pc in range(self.gwas_pcs):
						dataframe['PC' + str(pc+1)] = np.array(covs)[:,int(pc)]
					model = smf.ols('phenotype ~ genotype + ' + ' + '.join(['PC' + str(pc+1) for pc in range(self.gwas_pcs)]), data=dataframe).fit()
					t = float(model.params[1])/model.bse[1]
					p = scipy.stats.norm.sf(abs(t))*2
					return [model.params[1], model.bse[1], t, p] 

			if int_naught:
				if covs is None:
					if boot:				
						dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
						model = smf.ols(formula='phenotype ~ genotype - 1',data = dataframe).fit()
						t = float(model.params[0])/model.bse[0]
						p = scipy.stats.norm.sf(abs(t))*2
						boots = 500
						samples = [np.random.choice([i for i in len(gens)],len(gens)) for x in range(1,501)]
						boots = [np.sum(phens[sample]*gens[sample])/np.sum(gens[sample]**2) for sample in samples]						
						return [model.params[0], np.std(boots), t, p]
						
					else:
						dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
						model = smf.ols(formula='phenotype ~ genotype - 1',data = dataframe).fit()
						t = float(model.params[0])/model.bse[0]
						p = scipy.stats.norm.sf(abs(t))*2
						return [model.params[0], model.bse[0], t, p] 

				else:
					covs = np.array(covs)[:,0]
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens, 'covariates':covs})
					model = smf.ols(formula='phenotype ~ genotype + covariates - 1',data = dataframe).fit()
					t = float(model.params[0])/model.bse[0]
					p = scipy.stats.norm.sf(abs(t))*2
					return [model.params[0], model.bse[0], t, p] 
	
	#This function can be made redundant depending on how I recode the above fxn
	def estimate_effects(self, phenotypes, genotypes,covs = None, int_naught = False, boot = False):

		effects = np.apply_along_axis(self.extract_lm_coeffs,0,genotypes, phens = phenotypes, covs = covs, int_naught = int_naught, boot = boot)
		return effects.T
	
	#generate ancestral allele frequences for n.loci loci under neutral sfs with Ne = n.anc.
	def gen_neut_sfs(self, nloci, nanc):
		#set initial probabilities for observing a certain number alleles given the population size
		#being set equal to nanc
		probs_uw = np.ones(nanc)/[float(i+1) for i in range(nanc)]
		#Convert allele counts to their corresponding probabilities
		probs = probs_uw/np.sum(probs_uw)
		#calculate the cumulative distribution function for each number of alleles in a population of size
		#equal to nanc
		cdf = np.cumsum(probs)
		#draw a random allele frequency for each locus
		percs = np.random.uniform(0.01,0.99,size=nloci)
		#find where the loci frequency falls on the cdf and scale it according to the allele
		cdfind = np.array([(np.sum(cdf <= perc) + 1)/float(nanc) for perc in percs])
		return cdfind 

	def drift_afs(self, afs, drift_param):
		#drift the ancestral allele frequencies via a truncated normal.
		#supply a vector of allele frequncies and a drift parameter dp, and one realization
		#of the post-drift allele frequencies is provided.

		#determine the amount of drift that has occurred in a population at each site
		#which is paramaterized by the standing allele frequency
		afs = np.array(afs).astype(float)
		drift_sd = np.sqrt(drift_param*afs*(1-afs))
		drift = np.random.normal(0., drift_sd)
		#change the allele frequencies from those drawn from the afs according to their drift
		drifted = afs + drift
		#set all alleles that have fixed or been removed and set them to 1 or 0, respectively
		drifted[drifted < 0] = 0
		drifted[drifted > 1] = 1
		return drifted
	
	def offspring_genetic_sampler(self, parent_genos):
		parents1 = parent_genos[:int(parent_genos.shape[0]/2.)]
		parents2 = parent_genos[int(parent_genos.shape[0]/2.):]

		p1_max = 0
		p2_max = 0
		allmax = 0
		while p1_max != 1 or p2_max != 1:
			# sample = 
			p1 = np.random.binomial(1,parents1)
			p2 = np.random.binomial(1,parents2)
			total = p1+p2
			p1_max = np.max(p1)
			p2_max = np.max(p2)

		return p1+p2
		
	def offspring_generator(self, gens_par1, gens_par2):
		#given two parent genotypes, draw an offspring genotype
		gens_par1 = gens_par1/2
		gens_par2 = gens_par2/2
		all_genos = np.vstack((gens_par1,gens_par2))

		offspring_alleles = np.apply_along_axis(self.offspring_genetic_sampler,0,all_genos)
		# print(offspring_alleles.shape)
		# print(offspring_alleles.max(axis=0))

		return offspring_alleles

	def generate_effects(self, sigma_mat, n_loci, n_trait_loci, maf_loci, neutral_architecture):

		#set the relative scale of the variance in the effects by the inverse of the variance in their MAF
		#this captures the relationship where, generally, rarer variants have larger variants
		if neutral_architecture:
			maf_variance_scaling = np.ones((self.n_trait_loci))
		else:
			maf_variance_scaling = (1/np.sqrt(2*maf_loci[:self.n_trait_loci]*(1.0-maf_loci[:self.n_trait_loci])))

		direct_effects_trait_loci = np.random.normal(0,maf_variance_scaling,self.n_trait_loci)
		effects_trait_loci = np.vstack((direct_effects_trait_loci[:self.n_trait_loci],np.zeros(self.n_trait_loci,))).T

		#rescale the variance of the direct effects to be the user specified value
		if self.indirect_effect_variance == 0:
			effects_trait_loci[:,0] = ((effects_trait_loci[:,0] - effects_trait_loci[:,0].mean()))/np.std(effects_trait_loci[:,0])*np.sqrt(self.direct_effect_variance)
		else:
			effects_trait_loci[:,0] = ((effects_trait_loci[:,0] - effects_trait_loci[:,0].mean()))/np.std(effects_trait_loci[:,0])*np.sqrt(self.direct_effect_variance)
			effects_trait_loci[:,1] = effects_trait_loci[:,1] - np.mean(effects_trait_loci[:,1])
			effects_trait_loci = (effects_trait_loci@np.linalg.inv(np.linalg.cholesky(np.cov(effects_trait_loci.T))))#@np.linalg.cholesky(sigma_mat)
		
		#Draw a bunch of zeros to represent the loci that do not contribute to trait variance
		effects_non_trait_loci = np.zeros((n_loci - n_trait_loci,2))
		all_effects = np.vstack((effects_trait_loci,effects_non_trait_loci))
		
		return all_effects

	def calculate_true_se(self,betas,genos,true_phenotypes):
		#Given an effect and the data that it was generated for, calculate the standard error
		return np.array([np.sqrt(np.sum((true_phenotypes-(betas[variant] * genos[:,variant]))**2)/(genos.shape[0]-2.))/np.sqrt(np.sum((genos[:,variant]-np.mean(genos[:,variant]))**2)) for variant in range(genos.shape[1])])

	def genetic_sampler(self, af, ncohort):
		sample_max = 0
		sample_min = 1
		while sample_max !=2 and sample_min!=0:
			sample = np.random.binomial(2,af,size = ncohort)
			sample_max = np.max(sample)
			sample_min = np.max(sample)
			#come back here and cosample parents so that they have a total of at least 3 alleles to
		return sample

	def genetic_data_generator(self, cohort_afs, n_loci, n_cohort, all_effects):
		#Generate the parental genotypes given the cohort allele frequency spectrum
		#there is an issue downstream with the standard GWAS if the entire sample space of alleles isn't sampled
		#doing some rejection sampling here to insure that
		gens_parent1 = np.apply_along_axis(self.genetic_sampler, 0, cohort_afs.reshape(1,n_loci), n_cohort)
		gens_parent2 = np.apply_along_axis(self.genetic_sampler, 0, cohort_afs.reshape(1,n_loci), n_cohort)

		#generate two siblings per parental duo
		gens_sib1 = self.offspring_generator(gens_parent1, gens_parent2)	
		gens_sib2 = self.offspring_generator(gens_parent1, gens_parent2)
		#compute the genetic component (direct and indirect) of offspring traits in 
		indirect_trait = np.dot(gens_parent1 + gens_parent2,all_effects[:,1])
		direct_trait_sib1 = np.dot(gens_sib1,all_effects[:,0])
		direct_trait_sib2 = np.dot(gens_sib2,all_effects[:,0])
		#combine the direct and indirect into the total genetic variance
		genetic_trait_sib1 = direct_trait_sib1 + indirect_trait
		# genetic_trait_sib1 = (genetic_trait_sib1-genetic_trait_sib1.mean())/genetic_trait_sib1.std()
		genetic_trait_sib2 = direct_trait_sib2 + indirect_trait
		# genetic_trait_sib2 = (genetic_trait_sib2-genetic_trait_sib2.mean())/genetic_trait_sib2.std()
		return gens_parent1, gens_parent2, gens_sib1, gens_sib2, indirect_trait, direct_trait_sib1, direct_trait_sib2, genetic_trait_sib1, genetic_trait_sib2

	def environmental_data_generator(self, env_variance, env_covariance, n_cohort, direct_trait_sib1, direct_trait_sib2,
		genetic_trait_sib1, genetic_trait_sib2):

		#construct the covariance matrix for environmental effects between sibling
		env_cov_mat = np.array([[1,env_covariance],[env_covariance,1]])

		# print(env_cov_mat)
		#generate first sibling environmental effects and scale by to have the specified
		#proportion of the total variance
		env_sibs = np.random.multivariate_normal(np.array([0,0]),env_cov_mat, size = n_cohort)
		env_sibs = (env_sibs@np.linalg.inv(np.linalg.cholesky(np.cov(env_sibs.T))))@np.linalg.cholesky(env_cov_mat)

		#draw phenotypic contributions and scale
		env_scaler = self.env_variance*np.var(direct_trait_sib1)
		
		if np.sqrt(env_variance) != 0:
			env_sib1 = (env_sibs[:,0]/np.std(env_sibs[:,0]))*np.sqrt(env_scaler)
			env_sib2 = (env_sibs[:,1]/np.std(env_sibs[:,1]))*np.sqrt(env_scaler)
		else:
			env_sib1 = env_sibs[:,0]*0
			env_sib2 = env_sibs[:,1]*0

		trait_sib1 = genetic_trait_sib1 + env_sib1
		trait_sib2 = genetic_trait_sib2 + env_sib2

		return env_sib1, env_sib2, trait_sib1, trait_sib2

	def runner(self, n_loci, n_standard_gwas, n_sib_gwas, n_target_cohort, n_trait_loci, direct_effect_variance, 
		indirect_effect_variance, direct_indirect_covariance, drift_param, env_variance, env_covariance,
		boot = False, nanc = 100, neutral_architecture = False):
		
		#Generate a covariance matrix for simulating direct and indirect effects
		sigma_mat_dir_indir = np.array([[direct_effect_variance,direct_indirect_covariance * np.sqrt(direct_effect_variance) * np.sqrt(indirect_effect_variance)],
										[direct_indirect_covariance * np.sqrt(direct_effect_variance) * np.sqrt(indirect_effect_variance), indirect_effect_variance]])

		#Generate the neutral, ancestral sfs and let it drift
		gwas_cohort_afs = self.gen_neut_sfs(n_loci, nanc)
		prs_cohort_afs = self.drift_afs(gwas_cohort_afs, drift_param)		
		all_effects = self.generate_effects(sigma_mat_dir_indir, n_loci, n_trait_loci, gwas_cohort_afs, neutral_architecture)
		self.all_effects = all_effects
		#Assign direct and dynastic (indirect) effect to the loci
		#For the sibling GWAS, draw genotypes for parents (using the pop0 afs)
		#then from each set of parents, draw two offspring pairs
		gens_parent1_sibgwas, gens_parent2_sibgwas, gens_sib1_sibgwas, gens_sib2_sibgwas, indirect_trait_sibgwas, direct_trait_sib1_sibgwas, direct_trait_sib2_sibgwas, genetic_trait_sib1_sibgwas, genetic_trait_sib2_sibgwas = self.genetic_data_generator(gwas_cohort_afs, n_loci, n_sib_gwas, all_effects)

		#repreat the process for the target population that the PRS/PCA decomp will be performed on
		gens_parent1_prs, gens_parent2_prs, gens_sib1_prs, gens_sib2_prs, indirect_trait_prs, direct_trait_sib1_prs, direct_trait_sib2_prs, genetic_trait_sib1_prs, genetic_trait_sib2_prs = self.genetic_data_generator(prs_cohort_afs, n_loci, n_target_cohort, all_effects)

		#Draw environmental noise. This is entirely unshared environment;
		#there is no shared env for siblings beyond the indirect genetic effects.
		#We set the variance of the unshared env. term to be ~equal to that of
		#the variance due to direct genetic effects.
		env_sib1_sibgwas, env_sib2_sibgwas, trait_sib1_sibgwas, trait_sib2_sibgwas = self.environmental_data_generator(env_variance, env_covariance, n_sib_gwas, direct_trait_sib1_sibgwas, direct_trait_sib2_sibgwas, genetic_trait_sib1_sibgwas, genetic_trait_sib2_sibgwas)
		env_sib1_prs, env_sib2_prs, trait_sib1_prs, trait_sib2_prs = self.environmental_data_generator(env_variance, env_covariance, n_target_cohort, direct_trait_sib1_prs, direct_trait_sib2_prs, genetic_trait_sib1_prs, genetic_trait_sib2_prs)
		#We have generated data for both the sibling gwas and the cohort of individuals that will be included in the target PRS
		#cohort. Now we will generate the larger (by default), standard GWAS cohort
		gens_parent1_standardgwas, gens_parent2_standardgwas, gens_sib1_standardgwas, gens_sib2_standardgwas, indirect_trait_standardgwas, direct_trait_sib1_standardgwas, direct_trait_sib2_standardgwas, genetic_trait_sib1_standardgwas, genetic_trait_sib2_standardgwas = self.genetic_data_generator(gwas_cohort_afs, n_loci, n_standard_gwas, all_effects)
		env_sib1_standardgwas, env_sib2_standardgwas, trait_sib1_standardgwas, trait_sib2_standardgwas = self.environmental_data_generator(env_variance, env_covariance, n_standard_gwas, direct_trait_sib1_standardgwas, direct_trait_sib2_standardgwas, genetic_trait_sib1_standardgwas, genetic_trait_sib2_standardgwas)
		#Run the PCA combining populations 0 and 1, but using only one sibling
		#from each sibship. This is the PCA for correcting the GWAS
		pca_genotypes = gens_sib1_standardgwas
		allgwas_cohort = pca_genotypes

		pca = PCA()
		pca.fit_transform(pca_genotypes)
		pc_eigenvecs_gwas = pca.components_.T	
		pc_eigenvals_gwas = pca.explained_variance_.T.reshape(self.nloci,)
		covariates = np.multiply(np.matmul(pca_genotypes,pc_eigenvecs_gwas),pc_eigenvals_gwas)[:,:self.gwas_pcs]

		# GWAS in combined population, using only one sib from each sibship
		gwas_combined = self.estimate_effects(phenotypes = trait_sib1_standardgwas,
			genotypes = pca_genotypes, covs = covariates, int_naught = self.int_naught, boot = self.boot)
		gwas_combined = pd.DataFrame(gwas_combined, columns = ['beta','se','t','p'])
		no2counts = gwas_combined.loc[gwas_combined.isnull().any(axis=1)].index.tolist()
		# print(gwas_combined.isna())
		# gwas_combined = loc[no2counts] = [1.,1.,1.,1.1]
		gwas_combined = gwas_combined.drop(no2counts)
		# print(np.var(gwas_combined['se']))
		#sibling gwas
		sibgwas_uncorr_diff = self.estimate_effects(phenotypes = trait_sib1_sibgwas-trait_sib2_sibgwas,
			genotypes = gens_sib1_sibgwas-gens_sib2_sibgwas, int_naught = True, boot = boot)
		sibgwas_uncorr_diff = pd.DataFrame(sibgwas_uncorr_diff, columns = ['beta','se','t','p'])
		# sibgwas_uncorr_diff.loc[no2counts] = [1.,1.,1.,1.1]
		sibgwas_uncorr_diff = sibgwas_uncorr_diff.drop(no2counts)

		#run the pca combining populations 0 and 1, now with prs set
		snps_to_keep=gwas_combined.index.tolist()
		prs_pca_genotypes = gens_sib1_prs[:,snps_to_keep]
		# self.nloci = gwas_combined.shape[0]
		pca = PCA()
		pca.fit_transform(prs_pca_genotypes)
		prs_eigenvecs = pca.components_.T	
		prs_eigenvals = pca.explained_variance_.T

		#calculate the PRS from the combined GWAS - should result in a 1xn array
		#need to fix this to access a vector of effects
		prs_gwas_combined = np.matmul(prs_pca_genotypes, np.array(gwas_combined['beta']).reshape(len(snps_to_keep),1))
		#calculate PRS from the sibling GWAS
		prs_sibling = np.matmul(prs_pca_genotypes, np.array(sibgwas_uncorr_diff['beta']).reshape(len(snps_to_keep),1))
		trait_prs = trait_sib1_prs
		
		self.gwas_combined = gwas_combined.astype(float)
		self.sibgwas_uncorr_diff = sibgwas_uncorr_diff.astype(float)
		self.true_standard_beta = all_effects.sum(axis = 1)[snps_to_keep]
		self.true_standard_ses = self.calculate_true_se(self.true_standard_beta,pca_genotypes[:,snps_to_keep],trait_sib1_standardgwas)**2
		self.true_sib_beta = all_effects[snps_to_keep,0]
		self.true_sib_ses = self.calculate_true_se(self.true_sib_beta,(gens_sib1_sibgwas-gens_sib2_sibgwas)[:,snps_to_keep],trait_sib1_sibgwas-trait_sib2_sibgwas)
		return gwas_combined, sibgwas_uncorr_diff, prs_eigenvecs, prs_eigenvals, prs_gwas_combined, prs_sibling, trait_prs, gwas_cohort_afs

	def outputs(self):
		return {'gwas_combined':self.gwas_combined, 'sibgwas_uncorr':self.sibgwas_uncorr_diff,
		'target_eigenvalues':self.pc_sol_prs_values, 'target_eigenvectors':self.pc_sol_prs_vectors,
		'prs_gwas_combo':self.prs_gwas_comb, 'sib_prs':self.prs_sib,'trait_vals_prs':self.trait_prs,
		'af':self.af, 'env_covariance':self.env_covariance, 'env_variance':self.env_variance, 
		'true_standard_beta':self.true_standard_beta, 'true_standard_ses':self.true_standard_ses,
		'true_sib_beta':self.true_sib_beta, 'true_sib_ses':self.true_sib_ses}










