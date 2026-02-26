import warnings
# from pandas.errors import SettingWithCopyWarning
# warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import sys
import os
import time, sys, traceback, argparse
from pysnptools.snpreader import Bed
from pysnptools.standardizer import Unit
from estimation import *
from datetime import datetime
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


__version__ = '1.0.1'
MASTHEAD = "\n*********************************************************************************\n"
MASTHEAD += "* PGSUS\n"
MASTHEAD += "* SAD variance decomposition for polygenic scores\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2026 Samuel Pattillo Smith, Olivia Smith, Doc Edge, and Arbel Harpak\n"
MASTHEAD += "* University of Texas and University of Southern California\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************************\n"

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
# FILE PATHS
parser.add_argument("--genetic-file", type=str, default="/.", dest = 'genetic_file', help="PLINK .bed file with genotypes for prediction sample")
parser.add_argument("--pop-gwas-file", type=str, default="/.", dest = 'popgwas', help="Population GWAS summary statistics file path")
parser.add_argument("--sib-gwas-file", type=str, default="/.", dest = 'sibgwas', help="Sibling GWAS summary statistics file path")
parser.add_argument("--pvalue", type=float, default = 1, dest = 'pval', help="Max p-value threshold for SNP inclusion. Default is p<1 (all SNPs included).")
parser.add_argument('--threshold-list', type = str, dest='threshold_list', default=None, help="Instead of single p-value, a comma-separated list of p-value thresholds and perform one analysis for each threshold.")
parser.add_argument("--out", type=str, default = './', dest = 'outpath', help="Output directory. Default location is current directory.")
parser.add_argument("--outfile-label", type=str, default="", dest = 'outlabel', help="Output file stem")
parser.add_argument('--anc-data', type=str, dest = 'anc_data', default = 'support_files/SNPalleles_1000Genomes_allsites.txt.gz', help="By default will use the 1000 Genomes file in support_files directory.")
parser.add_argument('--block-bounds', type=str, dest = 'block_bounds', default = 'support_files/Pickrell_breakpoints_EUR.bed', help="By default will use the breakpoints file in support_files directory.")

# COLUMN HEADER OPTIONS
parser.add_argument("--ascertainment-set", type=str, default = 'gwas', dest = 'ascertainment', help="Which GWAS summary statistics to threshold on when determining which SNPs to include for PGS. Options are \'sibs\' or \'gwas\'.")
parser.add_argument("--chrom", type=str, default = 'CHR', dest = 'CHR', help="Column header with chromosome labels")
parser.add_argument("--pos", type=str, default = 'POS', dest = 'POS', help="Column header with position labels")
parser.add_argument("--chrom-pos", type=str, default = 'chrom.pos', dest = 'chr_pos', help="Column header with SNP IDs in format \'CHR:POS\'")
parser.add_argument("--pop-effect", type=str, default = 'BETA', dest = 'POPBETA', help="Effect column header for population GWAS. Default is \'BETA\'.")
parser.add_argument("--pop-se", type=str, default = 'se', dest = 'POPSE', help="Standard error column header for population GWAS. Default is \'se\'.")
parser.add_argument("--pop-pval-col", type=str, default = 'P', dest = 'POP_P', help="P-value column header for population GWAS. Default is \'P\'.")
parser.add_argument("--sib-effect", type=str, default = 'BETA', dest = 'SIBBETA', help="Effect column header for sibling GWAS. Default is \'BETA\'.")
parser.add_argument("--sib-se", type=str, default = 'se', dest = 'SIBSE', help="Standard error column header for sibling GWAS. Default is \'se\'.")
parser.add_argument("--pval-col", type=str, default = 'P', dest = 'P', help="P-value column header for sibling GWAS. Default is \'P\'.")

# DECOMPOSITION PARAMETERS
parser.add_argument('--nboots', type=int, dest = 'nboots', default = 100, help="Number of bootstraps to perform in estimating the standard errror of the isotropic inflation factor using the Deming regression framework. Default is 100.")
parser.add_argument("--eigenvals", type=str, default = None, dest = 'eigenvalues', help="File containing the eigenvalues of the prediction sample genotype matrix saved as an \".npy\" formatted file. If this flag is absent then a numpy implementation of PCA will be performed on the PLINK binary file provided as --genetic-file.")
parser.add_argument("--eigenvecs", type=str, default = None, dest = 'eigenvecs', help="File containing the eigenvectors of the prediction sample genotype matrix saved as an \".npy\" formatted file. If this flag is absent then a numpy implementation of PCA will be performed on the PLINK binary file provided as --genetic-file.")
parser.add_argument('--permutation-test', default=False, action=argparse.BooleanOptionalAction, dest = 'block_perm', help="Perform PC-wise permutation test. Otherwise, only the isotropic inflation factor will be estimated.")
parser.add_argument('--perm-pcs', type=int, dest = 'pcs_to_test', default = 15, help="Number of top PCs to be tested for the PC-wise permutation procedure. Default is 15.")
parser.add_argument('--nperm', type = int, dest = 'nperm', default = 1000, help="Number of permutations to perform for each PC-wise decomposition in constructing the permutation based null. Default is 1,000.")
parser.add_argument('--aperm', default = False, action=argparse.BooleanOptionalAction, dest = 'aperm', help="Perform an additional permutation procedure where the signs of the effect sizes in each block are flipped. This is a less computationally intensive than the block permutation procedure where the number of permutations is absolute. Default is False.")
parser.add_argument('--aperm-alpha', type = float, dest = 'aperm_alpha', default = 0.05, help="Sign-flipping permutation procedure will be performed if --aperm is set to True. This flag sets the alpha threshold for significance in this procedure. Default is 0.05.")
parser.add_argument('--c', type = float, dest = 'c', default = 0.1, help="The constant c used in the block permutation procedure to determine sensitivity threshold of adaptive permutaiton procedure. Default is 0.1.")
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
ascertainment = args.ascertainment.lower()
if (ascertainment not in ["sibs", "gwas"]):
	raise ValueError("--ascertainment option must be either \'sibs\' or \'gwas\'")
eigenvalues = args.eigenvalues
eigenvecs = args.eigenvecs
outpath = args.outpath
outlabel = args.outlabel

nboots = args.nboots
block_perm = args.block_perm
pcs_to_test = args.pcs_to_test
nperm = args.nperm
nboots = args.nboots
threshold_list = args.threshold_list
if threshold_list != None: # If there's just a p-value and not a threshold_list, set list to be that p-value
	threshold_list = args.threshold_list.split(',')
else:
	threshold_list = [pval]

aperm = args.aperm
aperm_alpha = args.aperm_alpha
c = args.c

# Check for presence of all files before beginning work
file_list = [genetic_file, popgwas, sibgwas, anc_data, block_bounds, outpath] #, genetic_file]
# Check that both or neither of eigenvalues are provided
if ((eigenvalues != None) & (eigenvecs == None)) | ((eigenvalues == None) & (eigenvecs != None)):
	raise ValueError("Both --eigenvals and --eigenvecs files must be provided to use pre-computed PCA. Please provide both file paths or provide neither and PCA will be performed.")
# If the files are provided, add them to the list of files to check
if (eigenvalues != None) & (eigenvecs != None):
    file_list.append(eigenvalues)
    file_list.append(eigenvecs)
for f in file_list:
	if os.path.exists(f) == False:
		raise FileNotFoundError(f'{f} does not exist or could not be opened.')

if __name__ == '__main__':

	args = parser.parse_args()
	if args.outlabel == "":
		raise ValueError('--outfile-label is required.')

	print(MASTHEAD)
	print(f"SAD decomposition started at {datetime.now()}")

	print(f"[args] Population GWAS file: {popgwas}")
	print(f"[args] Sibling GWAS file: {sibgwas}")
	print(f"[args] Genotype file: {genetic_file}")
	# print(f"[args] Reference allele file: {args.anc_data}")
	print(f"[args] Block bounds file: {block_bounds}")
	if threshold_list != [pval]:
		print(f"[args] Evaluating at p-value thresholds: {threshold_list}")
	else:
		print(f"[args] P-value threshold: {pval}")
	print(f"[args] Ascertainment dataset: {ascertainment}")
	print(f"[args] Output directory: {outpath}")
	print(f"[args] Output file stem: {outlabel}")
	
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
	anc_data = pd.read_csv(anc_data,sep = '\t', compression='infer')
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
	threshold = threshold_list[0]
	sad = estimate_components(block_bounds, pc_genotypes, popgwas[POPBETA].astype(float), popgwas[POPSE].astype(float), \
		sibgwas[SIBBETA].astype(float), sibgwas[SIBSE].astype(float), block_snps, asc_ps.astype(float), \
		float(threshold_list[0]), outpath, outlabel, CHR, POS, aperm, aperm_alpha, c, pc_lower_bound=100, eigenvecs = eigenvecs, eigenvalues = eigenvalues, \
		boot_se = nboots, block_perm = block_perm, pcs_to_test = pcs_to_test, nperm = nperm)
	
	log.log('Done.\n')
	final = sad.outputs()
	eigenvecs = final['eigenvecs']
	eigenvalues = final['eigenvalues']
	
	log.log('Alpha estimate: ' + str(final['alpha']) + ' (' + str(final['alpha_se']) + ')\n')	
	log.log('nSNPs passing ascertainment filter: ' + str(final['nsnp'])+'\n')	
	
	if outlabel != '':
		newfile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.alpha.txt','w')
		alphasefile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.alpha.se.txt','w')
		nsnpfile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.nsnp.txt','w')
		statPath = (f"{outpath}/{outlabel}.pval.{str(threshold)}.stats.txt")
		final['var_totals'].to_csv(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.variance.totals.txt', sep = '\t')

	else: ### SEEMS UNNECESSARY, SET DEFAULT OUTLABEL TO PGSUS AND REMOVE. CHECK PATH AS PART OF FILE FINDING ###
		newfile = open(outpath + '/pval.' + str(threshold) + '.alpha.txt','w')
		nsnpfile = open(outpath + '/pval.' + str(threshold) + '.nsnp.txt','w')	
		statPath = (f"{outpath}/{outlabel}.pval.{str(threshold)}.stats.txt")
		final['var_totals'].to_csv(outpath + '/pval.' + str(threshold) + '.variance.totals.txt', sep = '\t')
	
	newfile.write(str(final['alpha']))
	alphasefile.write(str(final['alpha_se']))
	nsnpfile.write(str(final['nsnp']))
	pd.DataFrame({'out_label':[outlabel], 'n_snp':[final['nsnp']], 'alpha':[final['alpha']], 'alpha_se':[final['alpha_se']], 'p_value':[threshold]}).to_csv(statPath, header=True, index=False, sep="\t")
	
	if len(threshold_list) > 1:
		for threshold in threshold_list[1:]:
			sad = estimate_components(block_bounds, pc_genotypes, popgwas[POPBETA].astype(float), popgwas[POPSE].astype(float), \
				sibgwas[SIBBETA].astype(float), sibgwas[SIBSE].astype(float), block_snps, asc_ps.astype(float), \
				float(threshold), outpath, outlabel, CHR, POS, aperm, aperm_alpha, c, pc_lower_bound=100, eigenvecs = eigenvecs, eigenvalues = eigenvalues, \
				boot_se = nboots, block_perm = block_perm, pcs_to_test = pcs_to_test, nperm = nperm)
			log.log('Done.\n')
			final = sad.outputs()

			log.log('Alpha estimate: ' + str(final['alpha']) + ' (' + str(final['alpha_se']) + ')\n')	
			log.log('nSNPs passing ascertainment filter: ' + str(final['nsnp'])+'\n')	
			
			if outlabel != '':
				newfile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.alpha.txt','w')
				alphasefile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.alpha.se.txt','w')
				nsnpfile = open(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.nsnp.txt','w')
				statPath = (f"{outpath}/{outlabel}.pval.{str(threshold)}.stats.txt")
				final['var_totals'].to_csv(outpath + '/' + outlabel + '.pval.' + str(threshold) + '.variance.totals.txt', sep = '\t')

			else: ### SEEMS UNNECESSARY, SET DEFAULT OUTLABEL TO PGSUS AND REMOVE. CHECK PATH AS PART OF FILE FINDING ###
				newfile = open(outpath + '/pval.' + str(threshold) + '.alpha.txt','w')
				nsnpfile = open(outpath + '/pval.' + str(threshold) + '.nsnp.txt','w')	
				statPath = (f"{outpath}/{outlabel}.pval.{str(threshold)}.stats.txt")
				final['var_totals'].to_csv(outpath + '/pval.' + str(threshold) + '.variance.totals.txt', sep = '\t')
			
			newfile.write(str(final['alpha']))
			alphasefile.write(str(final['alpha_se']))
			nsnpfile.write(str(final['nsnp']))
			pd.DataFrame({'out_label':[outlabel], 'n_snp':[final['nsnp']], 'alpha':[final['alpha']], 'alpha_se':[final['alpha_se']], 'p_value':[threshold]}).to_csv(statPath, header=True, index=False, sep="\t")
	
	log.log('SAD decomposition complete.\n')	
	print(f"SAD decomposition finished at {datetime.now()}")
