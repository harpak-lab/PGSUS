import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import sys
import time, sys, traceback, argparse
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
#This script then reads in the datasets, lines up the markers, and estimates
#	variance components along each principal component. Also returns sampling variation
#	via bootstrapping over Pickrell's ~1700 blocks and significance testing via #sign-flipping in each block.
################################################################################################


__version__ = '1.0.0'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* SAD variance decomposition for polygenic scores\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2024 Samuel Pattillo Smith, Doc Edge, and Arbel Harpak\n"
MASTHEAD += "* University of Texas and University of Southern California\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

class Logger(object):

    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        self.log_fh.write(msg)

#read in necessary arguments and filepaths, save as appropriate variables
parser = argparse.ArgumentParser()
parser.add_argument("--genetic-file", type=str, default = './', dest = 'genetic_file')
parser.add_argument("--pop-gwas-file", type=str, default = './', dest = 'popgwas')
parser.add_argument("--sib-gwas-file", type=str, default = './', dest = 'sibgwas')
parser.add_argument("--pvalue", type=float, default = 1, dest = 'pval')
parser.add_argument("--ascertainment-set", type=str, default = 'gwas', dest = 'ascertainment')
parser.add_argument("--chrom", type=str, default = 'CHR', dest = 'CHR')
parser.add_argument("--pos", type=str, default = 'POS', dest = 'POS')
parser.add_argument("--chrom-pos", type=str, default = 'chrom.pos', dest = 'chr_pos')
parser.add_argument("--pop-effect", type=str, default = 'BETA', dest = 'POPBETA')
parser.add_argument("--pop-se", type=str, default = 'se', dest = 'POPSE')
parser.add_argument("--sib-effect", type=str, default = 'BETA', dest = 'SIBBETA')
parser.add_argument("--sib-se", type=str, default = 'se', dest = 'SIBSE')
parser.add_argument("--pval-col", type=str, default = 'P', dest = 'P')
parser.add_argument("--pop-pval-col", type=str, default = 'P', dest = 'POP_P')

parser.add_argument('--nboots', type=int, dest = 'nboots', default = 100)
parser.add_argument("--eigenvals", type=str, default = None, dest = 'eigenvalues')
parser.add_argument("--eigenvecs", type=str, default = None, dest = 'eigenvecs')

parser.add_argument('--permutation-test', default=False, action=argparse.BooleanOptionalAction, dest = 'block_perm')

parser.add_argument('--perm-pcs', type=int, dest = 'pcs_to_test', default = 100)
parser.add_argument('--nperm', type = int, dest = 'nperm', default = 1000)

parser.add_argument("--outfile-label", type=str, default = '', dest = 'outlabel')
parser.add_argument("--out", type=str, default = './', dest = 'outpath')
parser.add_argument('--anc-data',type=str, dest = 'anc_data',default = 'support_files/SNPalleles_1000Genomes_allsites.txt.gz')
parser.add_argument('--block-bounds',type=str, dest = 'block_bounds',default = 'support_files/Pickrell_breakpoints_EUR.bed')

args = parser.parse_args()
genetic_file = args.genetic_file
popgwas = args.popgwas
sibgwas = args.sibgwas
pval = args.pval
CHR = args.CHR
POS = args.POS
chr_pos = args.chr_pos

POPBETA = args.POPBETA
POPSE = args.POPSE
SIBBETA = args.SIBBETA
SIBSE = args.SIBSE
POP_P = args.POP_P
P = args.P
anc_data = args.anc_data
block_bounds = args.block_bounds
ascertainment = args.ascertainment
eigenvalues = args.eigenvalues
eigenvecs = args.eigenvecs
outpath = args.outpath
outlabel = args.outlabel

nboots = args.nboots
block_perm = args.block_perm
pcs_to_test = args.pcs_to_test
nperm = args.nperm
nboots = args.nboots

if __name__ == '__main__':

	args = parser.parse_args()
	if args.outpath is None:
		raise ValueError('--outpath is required.')

	log = Logger(args.outpath+'/pval.' + str(pval) + '.log')

	defaults = vars(parser.parse_args(''))
	opts = vars(args)
	non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
	header = MASTHEAD
	header += "Call: \n"
	header += './sad.py \\\n'
	options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
	header += '\n'.join(options).replace('True','').replace('False','')
	header = header[0:-1]+'\n'
	log.log(header)
	log.log('Beginning analysis at {T}\n'.format(T=time.ctime()))
	start_time = time.time()

	if eigenvalues == None:
		log.log('No eigenvalues or eigenvectors were provided. PCA will be performed.\n')
	else:
		log.log('Eigenvalues and eigenvectors provided. Skipping PCA.\n')

	log.log('Reading support files...')
	#load in the allele info
	anc_data = pd.read_csv(anc_data,sep = '\t', compression = 'gzip')	
	anc_data['chrom.pos'] = anc_data['SNP'].astype(str)
	log.log('Done.\n')
	#read in genetic data
	log.log('Reading population GWAS statistics...')
	#read in the summary statistics for populations and siblings
	if '.gz' in popgwas:
		popgwas = pd.read_csv(popgwas, sep = '\t', compression = 'gzip')
	else:
		popgwas = pd.read_csv(popgwas, sep = '\t')

	popgwas['SE'] = popgwas[POPSE].astype(float)
	log.log('Done.\n')

	log.log('Reading sibling GWAS statistics...')
	if '.gz' in sibgwas:
		sibgwas = pd.read_csv(sibgwas, sep = '\t', compression = 'gzip')
	else:
		sibgwas = pd.read_csv(sibgwas, sep = '\t')

	log.log('Done.\n')
	#subset the population gwas results to just the sites that are in the target genotype data
	if genetic_file != './': 
		log.log('Reading target sample genotypes...')
		genotypes = Bed(genetic_file, count_A1 = True)
		loci = ['%g'%(y[0]) + ':' + str(y[2])[:-2] for y in genotypes.pos]

		popgwas = popgwas[popgwas[chr_pos].isin(loci)]
		sibgwas = sibgwas[sibgwas[chr_pos].isin(loci)].drop_duplicates(keep='first')
		overlap_snps = sibgwas.merge(popgwas, on = chr_pos,how = 'inner')[[chr_pos]]
		popgwas = popgwas[popgwas[chr_pos].isin(overlap_snps[chr_pos].tolist())].drop_duplicates(subset = chr_pos)
		sibgwas = sibgwas[sibgwas[chr_pos].isin(overlap_snps[chr_pos].tolist())].drop_duplicates(subset = chr_pos)
		statsnps = popgwas[chr_pos].tolist()
		statsnps = genotypes.sid_to_index(statsnps)
		
		pc_genotypes = genotypes.read().val[:,statsnps]	
		snps_nans = np.unique(np.argwhere(np.isnan(pc_genotypes))[:,1])
		pc_genotypes = np.delete(pc_genotypes, snps_nans, 1)
		statsnps = np.delete(statsnps,snps_nans,0)
		popgwas = popgwas.drop(snps_nans)
		popgwas = popgwas.reset_index(drop = True)
		sibgwas = sibgwas.drop(snps_nans)
		sibgwas = sibgwas.reset_index(drop = True)
	else:
		pc_genotypes = ''
	block_snps = popgwas[[CHR, POS, POP_P]]

	log.log('Done.\n')
	
	if ascertainment == 'sibs':
		asc_ps = sibgwas[P]
		stat_to_geno_df = sibgwas
	elif ascertainment == 'gwas':
		asc_ps = popgwas[POP_P]
		stat_to_geno_df = popgwas

	log.log('Beginning SAD decomposition...')
	sad = estimate_components(block_bounds, pc_genotypes, popgwas[POPBETA].astype(float), popgwas[POPSE].astype(float), \
		sibgwas[SIBBETA].astype(float), sibgwas[SIBSE].astype(float), block_snps, asc_ps.astype(float), \
		pval, outpath, outlabel, CHR, POS, pc_lower_bound=100, eigenvecs= eigenvecs, eigenvalues = eigenvalues, \
		boot_se = nboots, block_perm = block_perm, pcs_to_test = pcs_to_test, nperm = nperm)
	log.log('Done.\n')

	final = sad.outputs()
	
	log.log('Alpha estimate: ' + str(final['alpha']) + ' (' + str(final['alpha_se']) + ')\n')	
	log.log('nSNPs passing ascertainment filter: ' + str(final['nsnp'])+'\n')	
	
	if outlabel != '':
		newfile = open(outpath + '/' + outlabel + '.pval.' + str(pval) + '.alpha.txt','w')
		alphasefile = open(outpath + '/' + outlabel + '.pval.' + str(pval) + '.alpha.se.txt','w')
		nsnpfile = open(outpath + '/' + outlabel + '.pval.' + str(pval) + '.nsnp.txt','w')
		final['var_totals'].to_csv(outpath + '/' + outlabel + '.pval.' + str(pval) + '.variance.totals.txt', sep = '\t')
	else:
		newfile = open(outpath + '/pval.' + str(pval) + '.alpha.txt','w')
		nsnpfile = open(outpath + '/pval.' + str(pval) + '.nsnp.txt','w')	
		final['var_totals'].to_csv(outpath + '/pval.' + str(pval) + '.variance.totals.txt', sep = '\t')
	
	newfile.write(str(final['alpha']))
	alphasefile.write(str(final['alpha_se']))
	nsnpfile.write(str(final['nsnp']))
	log.log('SAD decomposition complete.\n')	









