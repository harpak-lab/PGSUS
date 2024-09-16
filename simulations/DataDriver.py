import numpy as np 
import pandas as pd 
import sys
import argparse
from generate_gwas_sumstats import GWAS
from estimate_alpha import SAD

import statsmodels.formula.api as sm
import pickle
from scipy.stats import sem
import os
import warnings
import glob

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

#Goal: Simulate unlinked loci that contribute to a trait, via 
#correlated direct and indirect effects. We simulate members of two populations
#simulated by a small amount of drift (optionally zero). A parental generation is created, and
#then parental pairs are formed who each have two children.
#In the offspring generation, we conduct a sib GWAS by comparing phenotype and genotype
#within sibships, and a standard GWAS by drawing one member of each sibship
#and adjusting for a specified number of PCs.
#In a second sample, we compute polygenic scores, compute PCs, and project 
#sib and standard polygenic scores on the PCs.

###################################################################
#Run a simulation
#Necessary flags detailed below
###################################################################
parser = argparse.ArgumentParser()
#file labeling
parser.add_argument("--out-label", type=str, default = 'out_test', dest = 'outlabel')
parser.add_argument("--sim-directory", type=str, default = None, dest = 'sim_dir')

#sample sizes and number of variants
parser.add_argument("--nsim", type=int, default = 1, dest = 'nsim')
parser.add_argument("--ngwaspcs", type=int, default = 1, dest = 'ngwaspcs')
parser.add_argument("--nanc", type=int, default = 100, dest = 'nanc')
parser.add_argument("--nloci", type=int, default = 400, dest = 'nloci')
parser.add_argument("--n-trait-loci", type=int, default = 200, dest = 'n_trait_loci')

#parameters for family data for the sibling GWAS and PRS target samples
parser.add_argument("--n-standard-gwas", type=int, default = 20000, dest = 'n_standard_gwas')
parser.add_argument("--n-sib-gwas", type=int, default = 10000, dest = 'n_sib_gwas')
parser.add_argument("--n-target-cohort", type=int, default = 1000, dest = 'n_target_cohort')

#parameters for variance components
parser.add_argument("--direct-variance", type=float, default = 1.0, dest = 'direct_effect_variance')
parser.add_argument("--indirect-variance", type=float, default = 0.0, dest = 'indirect_effect_variance')
parser.add_argument("--direct-indirect-covar", type=float, default = 0.0, dest = 'direct_indirect_covariance')

#parameters to toggle to understand their effect on the estimator, alpha
parser.add_argument("--drift-param", type=float, default = 0, dest = 'drift_param')
parser.add_argument("--env-variance", type=float, default = 0, dest = 'env_variance')
parser.add_argument("--env-covariance", type=float, default = 0, dest = 'env_covariance')
parser.add_argument("--ascp", type=str, default = '1.0', dest = 'ascp')
parser.add_argument("--pc-lower-bound", type=int, default = 100, dest = 'pc_lower_bound')
parser.add_argument('--neutral-architecture', default=False, action=argparse.BooleanOptionalAction, dest = 'neutral_architecture')
parser.add_argument('--component-estimate-error', default=False, action=argparse.BooleanOptionalAction, dest = 'comp_estimate_error')
parser.add_argument('--sib-se-multiplier', default=1., type=float, dest = 'sib_se_multiplier')

args = parser.parse_args()

#begin simulations
if not os.path.isdir('simulation_cache/'):
	os.system('mkdir simulation_cache')

if args.sim_dir:
	files = glob.glob('simulation_cache/' + args.sim_dir + '/sim.*.gwas.pkl')
	counter = 1
	for file in files:
		with open(file, 'rb') as f:
			simulation_gwas = pickle.load(f)
		decomposition = SAD(simulation_gwas, float(args.ascp), pc_lower_bound = args.pc_lower_bound)
		saddf = decomposition.outputs()
		outfile = file.replace('gwas','sad')
		with open(outfile, 'wb') as handle:
			pickle.dump(saddf, handle, protocol=pickle.HIGHEST_PROTOCOL)
		counter+=1

else:
	params = {'nsims':args.nsim, 'ngwaspcs':args.ngwaspcs, 'nanc':args.nanc, 'nloci':args.nloci, 'n_trait_loci':args.n_trait_loci, 
			  'n_standard_gwas':args.n_standard_gwas, 'n_sib_gwas':args.n_sib_gwas, 'n_target_cohort':args.n_target_cohort, 
			  'direct_effect_variance':args.direct_effect_variance, 'indirect_effect_variance':args.indirect_effect_variance, 
			  'direct_indirect_covariance':args.direct_indirect_covariance, 'drift_param':args.drift_param, 'env_variance':args.env_variance,
			  'env_covariance':args.env_covariance, 'ascp':args.ascp, 'pc_lower_bound':args.pc_lower_bound, 
			  'neutral_architecture':args.neutral_architecture}
	

	if args.neutral_architecture:
		architecture = 'neutral'
	else:
		architecture = 'maf'
	outlabel = 'sim.{}.architecture.pval.{}.npcs.{}.env.var.{}.env.cov.{}.nsib.{}.nstandard.{}.ntraitloci.{}'.format(architecture, args.ascp, args.ngwaspcs, args.env_variance, args.env_covariance, args.n_sib_gwas, args.n_standard_gwas, args.n_trait_loci)
		
	sim_params_manifest = pd.DataFrame.from_dict(params, orient = 'index')
	sim_params_manifest.to_csv('simulation_cache/sim.params.' + outlabel + '.txt',sep = '\t', index = False)

	print("Simulations for the following parameters:\n\
		N simluations: {}\n\
		N PCs in standard GWAS: {}\n\
		N ancestral population: {}\n\
		N total loci = {} \n\
		N trait loci = {} \n\
		N standard GWAS = {} \n\
		N sibling GWAS = {} \n\
		N target cohort = {} \n\
		Direct effect variance = {} \n\
		Indirect effect variance = {} \n\
		Direct-Indirect covariance = {} \n\
		Drift paramater = {} \n\
		Environmental variance = {} \n\
		Environmental covariance between siblings = {} \n\
		Ascertainment p-values = {} \n\
		Alpha PC lower bound = {} \n\
		Neutral architecture = {} \n \
		".format(args.nsim, args.ngwaspcs, args.nanc, args.nloci, args.n_trait_loci, args.n_standard_gwas, args.n_sib_gwas, 
			args.n_target_cohort, args.direct_effect_variance, args.indirect_effect_variance, args.direct_indirect_covariance, 
			args.drift_param, args.env_variance, args.env_covariance, args.ascp, args.pc_lower_bound, args.neutral_architecture))

	if os.path.exists('simulation_cache/' + outlabel):
		pass 
	else:
		os.system('mkdir simulation_cache/' + outlabel)

	for sim in range(1,args.nsim+1):
		simlabel = ''.join(["{}".format(np.random.randint(0, 9)) for num in range(0, 10)])
		simulation_gwas = GWAS(nloci=args.nloci, nanc=args.nanc, n_standard_gwas = args.n_standard_gwas, n_sib_gwas = args.n_sib_gwas, n_target_cohort = args.n_target_cohort,
				n_trait_loci = args.n_trait_loci, direct_effect_variance = args.direct_effect_variance, indirect_effect_variance = args.indirect_effect_variance,
				direct_indirect_covariance = args.direct_indirect_covariance, drift_param = args.drift_param, env_variance = args.env_variance, env_covariance = args.env_covariance,
				gwas_pcs = args.ngwaspcs, boot = False, neutral_architecture = args.neutral_architecture, covs = None)
		
		simulation_gwas = simulation_gwas.outputs()
		
		with open('simulation_cache/' + outlabel + '/sim.' + simlabel + '.gwas.pkl', 'wb') as handle:
			pickle.dump(simulation_gwas, handle, protocol=pickle.HIGHEST_PROTOCOL)

		print('Estimating alpha, Simulation: ' + str(sim))
		
		decomposition = SAD(simulation_gwas, float(args.ascp), pc_lower_bound = args.pc_lower_bound, sib_se_multiplier = args.sib_se_multiplier)
		saddf = decomposition.outputs()
		print(saddf['alpha'])

		with open('simulation_cache/' + outlabel + '/sim.' + simlabel + '.sad.pkl', 'wb') as handle:
			pickle.dump(saddf, handle, protocol=pickle.HIGHEST_PROTOCOL)



