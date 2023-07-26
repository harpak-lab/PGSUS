import pandas as pd
import numpy as np
import sys
import argparse
from pysnptools.snpreader import Bed
from pysnptools.standardizer import Unit
from estimation import *
################################################################################################
#This script is used to perform SAD variance decomposition in polygenic risk scores.
#The following inputs are necessary in order to perform SAD analysis:
#	- A genetic dataset for the target cohort (the prediction cohort where PCA is performed)
#	- A file containing summary statistics from a population level GWAS
#	- A file containing summary statistics form a sibling level GWAS
#	- A p-value threshold to determine which variants are used in the SAD decomposition
#	- The set of effect sizes in which to perform the ascertainment ('gwas' and 'sibs')
#	- A path to the directory where out files should be written
#
#IMPORTANT: the datasets with effect sizes (gwas and sib) need to be pre-processed to meet a few criteria. 
#	In particular, a column called "POS" that identifies loci by chromosome and position in GRCh37 coordinates.
#	All  effect sizes need to be in the same units (e.g. all in cm for a length measurement, or all in SD units)
#	The same  data on the same, in the same order, for a set of SNPs
#	Standard error of effect-size estimate named 'SE'
#	P-value for marker named 'PVAL'
#	Effect sizes polarized consistently across datasets (i.e. same alternate allele)
#	Effect-size estimate named 'BETA'
#
#For height, there is a script decomp_PRS_preprocess_height.R that performs these tasks; these can be 
#	performed prior to clumping. Clumping has its effect via the set of genotypes read in, 
#	and could be peformed either as part of this script or could precede it. 
#
#This script then reads in the datasets, lines up the markers, and estimates
#	variance components along each principal component. Also returns sampling variation
#	via bootstrapping over Pickrell's ~1700 blocks and significance testing via #sign-flipping in each block.
################################################################################################

#read in necessary arguments and filepaths, save as appropriate variables
parser = argparse.ArgumentParser()
parser.add_argument("--genetic-file", type=str, default = './', dest = 'genetic_file')
parser.add_argument("--pop-gwas-file", type=str, default = './', dest = 'popgwas')
parser.add_argument("--sib-gwas-file", type=str, default = './', dest = 'sibgwas')
parser.add_argument("--pvalue", type=float, default = 1, dest = 'pval')
parser.add_argument("--ascertainment-set", type=str, default = 'gwas', dest = 'ascertainment')
parser.add_argument("--nperm", type=int, default = 10000, dest = 'nperm')
parser.add_argument("--chrom", type=str, default = 'CHR', dest = 'CHR')
parser.add_argument("--pos", type=str, default = 'POS', dest = 'POS')
parser.add_argument("--pop-effect", type=str, default = 'BETA', dest = 'POPBETA')
parser.add_argument("--pop-stat", type=str, default = 'T_STAT', dest = 'POPSTAT')
parser.add_argument("--pval-col", type=str, default = 'P', dest = 'POP_P')
parser.add_argument("--out", type=str, default = './', dest = 'outpath')

args = parser.parse_args()
genetic_file = args.genetic_file
popgwas = args.popgwas
sibgwas = args.sibgwas
pval = args.pval
CHR = args.CHR
POS = args.POS
POPBETA = args.POPBETA
POPSTAT = args.POPSTAT
POP_P = args.POP_P
ascertainment = args.ascertainment
outpath = args.outpath

#load in the allele info
anc_data = pd.read_csv('support_files/SNPalleles_1000Genomes_allsitesMAF01.txt',sep = '\t')
anc_data['chrom.pos'] = anc_data['CHR'].astype(str) + '_' + anc_data['POS'].astype(str)

#read in genetic data

#read in the summary statistics for populations and siblings
popgwas = pd.read_csv(popgwas, sep = '\t', compression = 'gzip')
popgwas['SE'] = popgwas[POPBETA].astype(float)/popgwas[POPSTAT].astype(float)

sibgwas = pd.read_csv(sibgwas, sep = '\t', compression = 'gzip')

#subset the population gwas results to just the sites that are in the target genotype data
genotypes = Bed(genetic_file, count_A1 = True)
loci = ['%g'%(y[0]) + ':' + str(y[2])[:-2] for y in genotypes.pos]

popgwas = popgwas[popgwas['chrom.pos'].isin(loci)]
sibgwas = sibgwas[sibgwas['chrom.pos'].isin(loci)]

statsnps = popgwas['chrom.pos'].tolist()
statsnps = genotypes.sid_to_index(statsnps)

# pc_genotypes = Unit().standardize(genotypes.read().val)
pc_genotypes = genotypes.read().val

if ascertainment == 'sibs':
	asc_ps = sibgwas['EMP1']

elif ascertainment == 'gwas':
	asc_ps = popgwas[POP_P]


sad = estimate_components(pc_genotypes, popgwas['BETA'], popgwas['se'], sibgwas['BETA'], sibgwas['se'], asc_ps, pval, 
	pc_lower_bound=100, nperm = 10000, weight = True, eigenvecs= None, eigenvalues = None)

final = sad.outputs()
newfile = open(outpath + '/pval.' + str(pval) + '.gamma.txt','w')
sefile = open(outpath + '/pval.' + str(pval) + '.se.txt','w')
print(final['gamma'])
print(final['se'])
newfile.write(str(final['gamma']))
sefile.write(str(final['se']))










