import numpy as np
import pandas as pd
import sys
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()
from sklearn.decomposition import PCA
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

mass = importr('MASS')
base = importr('base')
stats = importr('stats')

r = robjects.r
FloatVector = robjects.FloatVector

class GWAS(object):

	def __init__(self, nloci, nanc, nfams_pop0, nfams_pop1, n_trait_loci, direct_effect_variance,
		indirect_effect_variance, direct_indirect_correlation ,drift_param_pop0, drift_param_pop1,
		env_diff, nfams_pop0_prs, nfams_pop1_prs, gwas_pcs, correct = True, boot = False, 
		scale_af = False, covs = None, int_naught = False):
		
		np.random.seed(8675309)
		print("GWAS performed with the following parameters:\n\
        	N total loci = {} \n\
        	N trait loci = {} \n\
        	N anc = {} \n\
        	N fams pop0 = {} \n\
        	N fams pop1 = {} \n\
        	N fams pop0 prs = {} \n\
        	N fams pop1 prs = {} \n\
        	Drift paramater pop0 = {} \n\
        	Drift paramater pop1 = {} \n\
        	Direct effect variance = {} \n\
        	Indirect effect variance = {} \n\
        	Direct-Indirect correlation = {} \n\
        	Environmental difference = {} \n\
        	Covariates = {} \n\
        	Bootstraps = {} \n\
        	Scale allele freq = {} \n \
        	".format(nloci, n_trait_loci, nanc, nfams_pop0, nfams_pop1, 
        							nfams_pop0_prs, nfams_pop1_prs, drift_param_pop0, drift_param_pop1, 
        							direct_effect_variance, indirect_effect_variance, direct_indirect_correlation, 
        							env_diff, covs, boot, scale_af))

		self.int_naught = int_naught
		self.boot = boot
		self.nloci = nloci
		self.gwas_pcs = gwas_pcs
		self.gwas_combined, self.sibgwas_uncorr_diff, self.pc_sol_prs_vectors, self.pc_sol_prs_values, self.prs_gwas_comb, self.prs_sib,self.trait_prs,self.af = self.generate_effects_and_pcs(nloci, nfams_pop0, nfams_pop1, n_trait_loci, direct_effect_variance, 
		indirect_effect_variance, direct_indirect_correlation, drift_param_pop0,drift_param_pop1,env_diff,
		nfams_pop0_prs = nfams_pop0_prs, nfams_pop1_prs = nfams_pop1_prs, correct=correct, boot = boot, nanc = nanc,
		scale_af = scale_af)
	
		self.outputs()

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
		percs = np.random.uniform(0.01,0.99,nloci)
		#find where the loci frequency falls on the cdf and scale it according to the allele
		#frequency cdf
		cdfind = [(np.sum(cdf <= perc) + 1)/float(nanc) for perc in percs]
		return cdfind 

	#drift the ancestral allele frequencies via a truncated normal.
	#supply a vector of allele frequncies and a drift parameter dp, and one realization
	#of the post-drift allele frequencies is provided.
	def drift_afs(self, afs, drift_param):
		#determine the amount of drift that has occurred in a population at each site
		#which is paramaterized by the standing allele frequency
		afs = np.array(afs).astype(float)
		drift_sd = np.sqrt(drift_param*afs*(1-afs))
		drift = np.array(r.rnorm(float(len(afs)),float(0), FloatVector(drift_sd)))
		#change the allele frequencies from those drawn from the afs according to their drift
		drifted = afs + drift

		#set all alleles that have fixed or been removed and set them to 1 or 0, respectively
		drifted[drifted < 0] = 0
		drifted[drifted > 1] = 1
		return drifted

	#given two parent genotypes, draw an offspring genotype
	def offspring_generator(self, gens_par1, gens_par2):
		gens_par1 = gens_par1/2		
		allele1 = np.array(np.random.binomial(1,gens_par1.flatten(),gens_par1.size).reshape(gens_par1.shape[0],gens_par1.shape[1]))
		gens_par2 = gens_par2/2
		allele2 = np.array(np.random.binomial(1,gens_par2.flatten(),gens_par2.size).reshape(gens_par2.shape[0],gens_par2.shape[1]))
		offspring_alleles = allele1 + allele2

		return offspring_alleles

	#Run a GWAS for one locus.
	#Gens and phens are vectors of genotypes and phenotypes.
	#optionally, the user can add covariates to the regression,
	#set an intercept of 0 (for sib regression), or get SEs via bootstrap
	#Bootstrap is done w/ 500 bootstrap samples, and the SE is *not* used for
	#calculating the p value at this moment, just for the SE estimate.
	def extract_lm_coeffs(self, gens, phens, covs = None, int_naught = False, boot = False):

		#check to see if SNP is monomorphic
		if np.var(gens) == 0:
			 return np.array([0,np.nan,np.nan,1])
		else:

			#otherwise run regression, either without covariates or with them
			output_sumstats = []
			if int_naught == False:
				if covs is None:
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
					pandas2ri.activate()
					robjects.globalenv['dataframe'] = dataframe
					M = stats.lm('phenotype ~ genotype', data=base.as_symbol('dataframe'))
					lmsummary = base.summary(M).rx2('coefficients')
					return lmsummary[1,:]

				else:
					# covs = np.array(covs)[:,0]
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
					for pc in range(self.gwas_pcs):
						dataframe['PC' + str(pc+1)] = np.array(covs)[:,int(pc)]

					pandas2ri.activate()
					robjects.globalenv['dataframe'] = dataframe

					M = stats.lm('phenotype ~ genotype + ' + '+ '.join(['PC' + str(pc+1) for pc in range(self.gwas_pcs)]), data=base.as_symbol('dataframe'))

					lmsummary = base.summary(M).rx2('coefficients')
					return lmsummary[1,:]


			if int_naught:
				if covs is None:
					if boot:				

						dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
						pandas2ri.activate()
						robjects.globalenv['dataframe'] = dataframe
						M = stats.lm(formula='phenotype ~ genotype - 1',data = base.as_symbol('dataframe'))
						lmsummary = base.summary(M).rx2('coefficients')[0,:]
						boots = 500
						# lmsummary[1,:] = 
						samples = [np.random.choice([i for i in len(gens)],len(gens)) for x in range(1,501)]
						boots = [np.sum(phens[sample]*gens[sample])/np.sum(gens[sample]**2) for sample in samples] 
						#need to reset the standard error here using the follwoing
						lmsummary[1] = np.std(boots)
						return lmsummary
						#output regression with intercept, no covariates, no bootstraps
					else:

						dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens})
						pandas2ri.activate()
						robjects.globalenv['dataframe'] = dataframe
						M = stats.lm(formula='phenotype ~ genotype - 1',data = base.as_symbol('dataframe'))
						lmsummary = base.summary(M).rx2('coefficients')
						return lmsummary[0,:]

				else:
					covs = np.array(covs)[:,0]
					dataframe = pd.DataFrame({'genotype': gens, 'phenotype': phens, 'covariates':covs})
					pandas2ri.activate()
					robjects.globalenv['dataframe'] = dataframe

					M = stats.lm(formula='phenotype ~ genotype + covariates - 1',data = base.as_symbol('dataframe'))
					lmsummary = base.summary(M).rx2('coefficients')
					return lmsummary[0,:]
	
	#This function can be made redundant depending on how I recode the above fxn
	def estimate_effects(self,phenotypes, genotypes,covs = None, int_naught = False, boot = False):
		effects = np.apply_along_axis(self.extract_lm_coeffs,0, genotypes, phens = phenotypes, covs = covs, int_naught = int_naught, boot = boot)
		
		return effects.T

	def generate_effects_and_pcs(self, nloci, nfams_pop0, nfams_pop1, n_trait_loci, direct_effect_variance, 
		indirect_effect_variance, direct_indirect_correlation, drift_param_pop0,drift_param_pop1,env_diff,
		nfams_pop0_prs, nfams_pop1_prs, correct=True, boot = False, nanc = 100,
		scale_af = False):
		
		#Generate a covariance matrix for simulating direct and indirect effects
		sigma_mat_dir_indir = np.array([[direct_effect_variance,direct_indirect_correlation * np.sqrt(direct_effect_variance) * np.sqrt(indirect_effect_variance)],
										[direct_indirect_correlation * np.sqrt(direct_effect_variance) * np.sqrt(indirect_effect_variance), indirect_effect_variance]])
		
		
		#Generate the neutral, ancestral sfs and let it drift
		ancestral_afs = self.gen_neut_sfs(nloci,nanc)
		pop0_afs = self.drift_afs(ancestral_afs, drift_param_pop0)
		pop1_afs = self.drift_afs(ancestral_afs, drift_param_pop1)
		r_covar_vals = robjects.r.matrix(FloatVector(sigma_mat_dir_indir.ravel()), nrow = sigma_mat_dir_indir.shape[0])

		#Assign direct and dynastic (indirect) effect to the loci
		effects_trait_loci = np.array(mass.mvrnorm(float(n_trait_loci),FloatVector(np.zeros(2)),r_covar_vals)).reshape(n_trait_loci,2)
		effects_non_trait_loci = np.zeros((nloci - n_trait_loci,2))
		effects_loci = np.vstack((effects_trait_loci,effects_non_trait_loci))

		#In population 0, draw genotypes from parents (using the pop0 afs)
		#then from each set of parents, draw two offspring pairs
		gens_parent1_pop0 = np.array(r.matrix(r.rbinom(float(nfams_pop0)*float(nloci),2.0, r.rep(FloatVector(pop0_afs),float(nfams_pop0))), ncol = nfams_pop0))
		gens_parent2_pop0 = np.array(r.matrix(r.rbinom(float(nfams_pop0)*float(nloci),2.0, r.rep(FloatVector(pop0_afs),float(nfams_pop0))), ncol = nfams_pop0))
		gens_sib1_pop0 = self.offspring_generator(gens_parent1_pop0, gens_parent2_pop0)	
		gens_sib2_pop0 = self.offspring_generator(gens_parent1_pop0, gens_parent2_pop0)

		#repreat the process for the population that the PRS/PCA decomp will be performed on
		gens_parent1_pop0_prs = np.array(r.matrix(r.rbinom(float(nfams_pop0_prs)*float(nloci),2.0, r.rep(FloatVector(pop0_afs),float(nfams_pop0_prs))), ncol = nfams_pop0_prs))
		gens_parent2_pop0_prs = np.array(r.matrix(r.rbinom(float(nfams_pop0_prs)*float(nloci),2.0, r.rep(FloatVector(pop0_afs),float(nfams_pop0_prs))), ncol = nfams_pop0_prs))
		gens_os_pop0_prs = self.offspring_generator(gens_parent1_pop0_prs, gens_parent2_pop0_prs)
				
		#compute the genetic component (direct and indirect) of offspring traits in 
		#pop 0 - let column 1 be the direct effects
		indirect_trait_pop0 = np.matmul((gens_parent1_pop0 + gens_parent2_pop0).T,effects_loci[:,1])
		direct_trait_pop0_sib1 = np.matmul(gens_sib1_pop0.T,effects_loci[:,0])
		direct_trait_pop0_sib2 = np.matmul(gens_sib2_pop0.T,effects_loci[:,0])
		#combine the direct and indirect into the total genetic variance
		genetic_trait_pop0_sib1 = direct_trait_pop0_sib1 + indirect_trait_pop0
		genetic_trait_pop0_sib2 = direct_trait_pop0_sib2 + indirect_trait_pop0

		#repeat the process for the fams that are in the PRS prediction cohort
		indirect_trait_pop0_prs = np.matmul((gens_parent1_pop0_prs+gens_parent2_pop0_prs).T,effects_loci[:,1])
		direct_trait_pop0_prs = np.matmul(gens_os_pop0_prs.T,effects_loci[:,0])
		genetic_trait_pop0_prs = indirect_trait_pop0_prs + direct_trait_pop0_prs

		#repeat the entire process for population 1
		gens_parent1_pop1 = np.array(r.matrix(r.rbinom(float(nfams_pop1)*float(nloci),2.0, r.rep(FloatVector(pop1_afs),float(nfams_pop1))), ncol = nfams_pop1))
		gens_parent2_pop1 = np.array(r.matrix(r.rbinom(float(nfams_pop1)*float(nloci),2.0, r.rep(FloatVector(pop1_afs),float(nfams_pop1))), ncol = nfams_pop1))
		gens_sib1_pop1 = self.offspring_generator(gens_parent1_pop1,gens_parent2_pop1)
		gens_sib2_pop1 = self.offspring_generator(gens_parent1_pop1,gens_parent2_pop1)
		
		#get the variance components in pop1
		gens_parent1_pop1_prs = np.array(r.matrix(r.rbinom(float(nfams_pop1_prs)*float(nloci),2.0, r.rep(FloatVector(pop1_afs),float(nfams_pop1_prs))), ncol = nfams_pop1_prs))
		gens_parent2_pop1_prs = np.array(r.matrix(r.rbinom(float(nfams_pop1_prs)*float(nloci),2.0, r.rep(FloatVector(pop1_afs),float(nfams_pop1_prs))), ncol = nfams_pop1_prs))
		gens_os_pop1_prs = self.offspring_generator(gens_parent1_pop1_prs,gens_parent2_pop1_prs)
		
		#compute the genetic component (direct and indirect) of offspring traits in 
		#pop 1 - let column 2 be the direct effects
		indirect_trait_pop1 = np.matmul((gens_parent1_pop1+gens_parent2_pop1).T,effects_loci[:,1])
		direct_trait_pop1_sib1 = np.matmul(gens_sib1_pop1.T,effects_loci[:,0])
		direct_trait_pop1_sib2 = np.matmul(gens_sib2_pop1.T,effects_loci[:,0])
		#combine for total genetic variance
		genetic_trait_pop1_sib1 = direct_trait_pop1_sib1 + indirect_trait_pop1
		genetic_trait_pop1_sib2 = direct_trait_pop1_sib2 + indirect_trait_pop1
		#repeat for pop1 - prs cohort
		#calculte prs for pop1 cohort
		indirect_trait_pop1_prs = np.matmul((gens_parent1_pop1_prs + gens_parent2_pop1_prs).T,effects_loci[:,1])
		direct_trait_pop1_prs = np.matmul(gens_os_pop1_prs.T,effects_loci[:,0])
		genetic_trait_pop1_prs = indirect_trait_pop1_prs + direct_trait_pop1_prs

		################################################################
		#Draw environmental noise. This is entirely unshared environment;
		#there is no shared env for siblings beyond the indirect genetic effects.
		#We set the variance of the unshared env. term to be ~equal to that of
		#the variance due to direct genetic effects.
		################################################################
		env_pop0_sib1 = np.array(r.rnorm(float(nfams_pop0),float(0),float(1)))
		if np.std(direct_trait_pop0_sib1) != 0:
			env_pop0_sib1 = env_pop0_sib1*np.std(direct_trait_pop0_sib1)

		env_pop0_sib2 = np.array(r.rnorm(float(nfams_pop0),float(0),float(1)))
		if np.std(direct_trait_pop0_sib2) != 0:
			env_pop0_sib2 = env_pop0_sib2*np.std(direct_trait_pop0_sib2)

		env_pop0_prs = np.array(r.rnorm(float(nfams_pop0),float(0),float(1)))
		if np.std(direct_trait_pop0_prs) != 0:
			env_pop0_prs = env_pop0_prs*np.std(direct_trait_pop0_prs)

		env_pop1_sib1 = np.array(r.rnorm(float(nfams_pop1),float(0),float(1)))
		if np.std(direct_trait_pop1_sib1) != 0:
			env_pop1_sib1 = env_pop1_sib1*np.std(direct_trait_pop1_sib1)

		env_pop1_sib2 = np.array(r.rnorm(float(nfams_pop1),float(0),float(1)))
		if np.std(direct_trait_pop1_sib2) != 0:
			env_pop1_sib2 = env_pop1_sib2*np.std(direct_trait_pop1_sib2)

		env_pop1_prs = np.array(r.rnorm(float(nfams_pop1),float(0),float(1)))
		if np.std(direct_trait_pop1_prs) != 0:
			env_pop1_prs = env_pop1_prs*np.std(direct_trait_pop1_prs)

		#Compute the trait value for each pair of sibings in each population
		trait_pop0_sib1 = genetic_trait_pop0_sib1 + env_pop0_sib1
		trait_pop0_sib2 = genetic_trait_pop0_sib2 + env_pop0_sib2
		trait_pop0_prs = genetic_trait_pop0_prs + env_pop0_prs

		trait_pop1_sib1 = genetic_trait_pop1_sib1 + env_pop1_sib1
		trait_pop1_sib2 = genetic_trait_pop1_sib2 + env_pop1_sib2
		trait_pop1_prs = genetic_trait_pop1_prs + env_pop1_prs

		#Run the PCA combining populations 0 and 1, but using only one sibling
		#from each sibship. This is the PCA for correcting the GWAS
		pca_genotypes = np.hstack((gens_sib1_pop0,gens_sib1_pop1)).T
		
		if scale_af is False:
			#need to check dimensions here
			allgwas_cohort = np.hstack((gens_sib1_pop0,gens_sib1_pop1)).T
			matnorm = r.matrix(FloatVector(allgwas_cohort.flatten(order='K')),nrow = allgwas_cohort.shape[0])
			pc_sol = r.prcomp(matnorm)
			pc_sol = pc_sol.rx2('x')

		elif scale_af is True:
			allgwas_cohort = np.hstack((gens_sib1_pop0,gens_sib1_pop1)).T
			af = np.mean(allgwas_cohort, axis = 0)/float(2)
			# af = af.fillna(0)
			# print(np.where(np.isinf(allgwas_cohort)))
			allgwas_cohort = (allgwas_cohort - 2*af)/np.sqrt(2*af*(1-af))
			matnorm = r.matrix(FloatVector(allgwas_cohort.flatten(order='K')),nrow = allgwas_cohort.shape[0])
			pc_sol = r.prcomp(matnorm)
			pc_sol = pc_sol.rx2('x')
			
		# GWAS in combined population, using only one sib from each sibship
		if correct:
			gwas_combined = self.estimate_effects(phenotypes = np.concatenate((trait_pop0_sib1,trait_pop1_sib1)),
				genotypes = pca_genotypes, covs = pc_sol, int_naught = self.int_naught, boot = self.boot)

		if not correct:
			gwas_combined = self.estimate_effects(phenotypes = np.concatenate((trait_pop0_sib1.T,trait_pop1_sib1.T)),
				genotypes = pca_genotypes, int_naught = self.int_naught, boot = self.boot)

		gwas_combined = pd.DataFrame(gwas_combined, columns = ['beta','se','t','p'])
		
		#sibling gwas
		sibgwas_uncorr_diff = self.estimate_effects(phenotypes = np.concatenate((trait_pop0_sib1-trait_pop0_sib2,trait_pop1_sib1-trait_pop1_sib2)),
			genotypes = np.hstack((gens_sib1_pop0-gens_sib2_pop0,gens_sib1_pop1-gens_sib2_pop1)).T, int_naught = True, boot = boot)

		sibgwas_uncorr_diff = pd.DataFrame(sibgwas_uncorr_diff, columns = ['beta','se','t','p'])

		pc_sol_gwas = pc_sol

		#run the pca combining populations 0 and 1, now with prs set
		prs_pca_genotypes = np.hstack((gens_os_pop0_prs,gens_os_pop1_prs)).T

		# af = np.mean(prs_pca_genotypes,axis =0)/2
		if scale_af is False:
			matnorm = r.matrix(FloatVector(prs_pca_genotypes.flatten(order='K')),nrow = prs_pca_genotypes.shape[0])
			pc_sol = r.prcomp(matnorm)
			pc_sol_prs_eigenvalues = np.array(pc_sol.rx2('sdev'))**2
			pc_sol_prs_eigenvectors = np.array(pc_sol.rx2('rotation'))
			pc_sol = pc_sol.rx2('x')

		elif scale_af is True:

			allgwas_cohort = np.hstack((gens_os_pop0_prs,gens_os_pop1_prs)).T
			af = np.mean(allgwas_cohort, axis = 0)/float(2)
			allgwas_cohort = (allgwas_cohort - 2*af)/np.sqrt(2*af*(1-af))

			matnorm = r.matrix(FloatVector(allgwas_cohort.flatten(order='K')),nrow = allgwas_cohort.shape[0])
			pc_sol = r.prcomp(matnorm)
			pc_sol_prs_eigenvalues = np.array(pc_sol.rx2('sdev'))**2
			pc_sol_prs_eigenvectors = np.array(pc_sol.rx2('rotation'))
			pc_sol = pc_sol.rx2('x')

		#calculate the PRS from the combined GWAS - should result in a 1xn array
		#need to fix this to access a vector of effects
		prs_gwas_combined = np.matmul(prs_pca_genotypes, np.array(gwas_combined['beta']).reshape(self.nloci,1))
		
		#calculate PRS from the sibling GWAS
		prs_sibling = np.matmul(prs_pca_genotypes, np.array(sibgwas_uncorr_diff['beta']).reshape(self.nloci,1))
		trait_prs = np.vstack((trait_pop0_prs,trait_pop1_prs))
		
		gwas_combined = gwas_combined.astype(float)
		sibgwas_uncorr_diff = sibgwas_uncorr_diff.astype(float)
		
		if scale_af:
			gwas_combined['af'] = af
			gwas_combined = gwas_combined.fillna(0)
			
			gwas_combined['beta'] = gwas_combined['beta']*np.sqrt(2*gwas_combined['af']*(1-gwas_combined['af']))
			gwas_combined['se'] = gwas_combined['se']*np.sqrt(2*af*(1-af))

			sibgwas_uncorr_diff['af'] = af
			sibgwas_uncorr_diff = sibgwas_uncorr_diff.fillna(0)
			
			sibgwas_uncorr_diff['beta'] = sibgwas_uncorr_diff['beta']*np.sqrt(2*sibgwas_uncorr_diff['af']*(1-sibgwas_uncorr_diff['af']))
			sibgwas_uncorr_diff['se'] = sibgwas_uncorr_diff['se']*np.sqrt(2*sibgwas_uncorr_diff['af']*(1-sibgwas_uncorr_diff['af']))

		return gwas_combined, sibgwas_uncorr_diff, pc_sol_prs_eigenvectors, pc_sol_prs_eigenvalues, prs_gwas_combined, prs_sibling, trait_prs, af

	def outputs(self):
		return {'gwas_combined':self.gwas_combined, 'sibgwas_uncorr':self.sibgwas_uncorr_diff,
		'prs_pcs_vectors':self.pc_sol_prs_vectors, 'prs_pcs_values':self.pc_sol_prs_values,
		'prs_gwas_combo':self.prs_gwas_comb, 'sib_prs':self.prs_sib,'trait_vals_prs':self.trait_prs,'af':self.af}










