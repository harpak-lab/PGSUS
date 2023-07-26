import numpy as np 
import pandas as pd 
import sys
import argparse
from generate_gwas_sumstats import GWAS
from estimate_gamma import SAD
import statsmodels.formula.api as sm
import pickle
from scipy.stats import sem
import os

#June 2023 adapted from April 2023, by Doc Edge
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
parser.add_argument("--nsim", type=int, default = 1, dest = 'nsim')
parser.add_argument("--ngwaspcs", type=int, default = 1, dest = 'ngwaspcs')
parser.add_argument("--nanc", type=int, default = 100, dest = 'nanc')
parser.add_argument("--nloci", type=int, default = 400, dest = 'nloci')
parser.add_argument("--ntraitloci", type=int, default = 200, dest = 'ntraitloci')
parser.add_argument("--nfams-pop0", type=int, default = 2000, dest = 'nfams_pop0')
parser.add_argument("--nfams-pop1", type=int, default = 2000, dest = 'nfams_pop1')
parser.add_argument("--nfams-pop0-prs", type=int, default = 2000, dest = 'nfams_pop0_prs')
parser.add_argument("--nfams-pop1-prs", type=int, default = 2000, dest = 'nfams_pop1_prs')
parser.add_argument("--direct-variance", type=float, default = 0.1, dest = 'direct_effect_variance')
parser.add_argument("--indirect-variance", type=float, default = 0.0, dest = 'indirect_effect_variance')
parser.add_argument("--direct-indirect-corr", type=float, default = 0.0, dest = 'direct_indirect_correlation')
parser.add_argument("--drift-pop0", type=float, default = 0, dest = 'drift_pop0')
parser.add_argument("--drift-pop1", type=float, default = 0, dest = 'drift_pop1')
parser.add_argument("--env-diff", type=float, default = 0, dest = 'env_difference')
parser.add_argument("--ascp", type=float, default = 1, dest = 'ascp')
parser.add_argument("--pc-lower-bound", type=int, default = 100, dest = 'pc_lower_bound')
parser.add_argument('--correct', default=False, action=argparse.BooleanOptionalAction, dest = 'correct')
parser.add_argument('--scale-af', default=True, action=argparse.BooleanOptionalAction, dest = 'scale_af')
args = parser.parse_args()


#begin simulations
gammas_basic_unweighted = np.zeros(args.nsim)
gammas_basic_weighted = np.zeros(args.nsim)

for sim in range(1,args.nsim+1):
	print('Simulation: ' + str(sim))
	if os.path.isfile('pythonsims/sim.' + str(sim) + '.' + str(args.ngwaspcs) + '.pc.pkl'):
		with open('pythonsims/sim.' + str(sim) + '.' + str(args.ngwaspcs) + '.pc.pkl', 'rb') as handle:
			simulation_gwas = pickle.load(handle)

	else:
		simulation_gwas = GWAS(nloci=args.nloci, nanc=args.nanc, nfams_pop0 = args.nfams_pop0, nfams_pop1 = args.nfams_pop1, n_trait_loci = args.ntraitloci,
				direct_effect_variance = args.direct_effect_variance, indirect_effect_variance = args.indirect_effect_variance,
				direct_indirect_correlation = args.direct_indirect_correlation, drift_param_pop0 = args.drift_pop0, drift_param_pop1 = args.drift_pop1, 
				env_diff = args.env_difference, nfams_pop0_prs = args.nfams_pop0_prs, nfams_pop1_prs = args.nfams_pop1_prs,
				gwas_pcs = args.ngwaspcs, correct=args.correct, boot = False, scale_af = args.scale_af, covs = None)
		simulation_gwas = simulation_gwas.outputs()

		with open('pythonsims/sim.' + str(sim) + '.' + str(args.ngwaspcs) + '.pc.pkl', 'wb') as handle:
			pickle.dump(simulation_gwas, handle, protocol=pickle.HIGHEST_PROTOCOL)
	
	decomposition = SAD(simulation_gwas, args.ascp, pc_lower_bound = args.pc_lower_bound)
	final = decomposition.outputs()
	gammas_basic_unweighted[sim-1] = final['gamma']

	decomposition = SAD(simulation_gwas, args.ascp, pc_lower_bound = args.pc_lower_bound, weight = True)
	final = decomposition.outputs()
	gammas_basic_weighted[sim-1] = final['gamma']

if not os.path.isdir('gamma_outputs/'):
	os.system('mkdir gamma_outputs')

pd.DataFrame(gammas_basic_unweighted.T).to_csv('gamma_outputs/unweighted.gammas.' + str(args.ngwaspcs) + '.gwaspcs.' + str(args.pc_lower_bound) + '.pc.lower.bound.' + str(args.ascp) + '.asc.p.txt',sep = '\t', index = False, header = False)
pd.DataFrame(gammas_basic_weighted.T).to_csv('gamma_outputs/weighted.gammas.' + str(args.ngwaspcs) + '.gwaspcs.' + str(args.pc_lower_bound) + '.pc.lower.bound.' + str(args.ascp) + '.asc.p.txt',sep = '\t', index = False, header = False)

print('Unweighted mean (python):' + str(np.mean(gammas_basic_unweighted)))
print('Unweighted std. err. (python):' + str(sem(gammas_basic_unweighted)))
print('Weighted mean (python):' + str(np.mean(gammas_basic_weighted)))
print('Weighted std. err. (python):' + str(sem(gammas_basic_unweighted)))


