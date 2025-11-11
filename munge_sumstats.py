import pandas as pd 
import numpy as np 
import time, sys, traceback, argparse
import os, subprocess
from scipy import stats
import warnings
from datetime import datetime
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

parser = argparse.ArgumentParser()
# parser.add_argument("--bfile", type=str, default = 'support_files/eur_1000G.noduplicates.maf01.snpsonly', dest = 'genetic_file')
parser.add_argument("--anc-data", type=str, default = 'support_files/SNPalleles_1000Genomes_allsites.txt.gz', dest = 'ancdata')
parser.add_argument("--pop-gwas-file", type=str, default = './', dest = 'popgwas')
parser.add_argument("--sib-perm-file", type=str, default = './', dest = 'sibgwasperm')
parser.add_argument("--chr", type=str, default = 'CHR', dest = 'chrom')
parser.add_argument("--pos", type=str, default = 'POS', dest = 'pos')
parser.add_argument("--snp-id", type=str, default = 'SNP', dest = 'snpid')
parser.add_argument('--outdir', type = str, dest = 'outdir')
parser.add_argument('--outlabel', type = str, dest = 'outlabel', default = '')
parser.add_argument('--standard-beta', type = str, dest = 'standard_beta', default = 'BETA')
parser.add_argument('--sib-beta', type = str, dest = 'sib_beta', default = 'BETA')
parser.add_argument('--pval', type = str, dest = 'P', default = 'P')
parser.add_argument('--preselected-snps', type=str, default = None, dest = 'snpset')
parser.add_argument('--a1', type = str, default = 'A1', dest='alt_allele')
parser.add_argument('--p-is-log', default=False, action=argparse.BooleanOptionalAction, dest = 'log10p')
parser.add_argument('--logp-col', default=None, dest = 'logp_col', type = str)

args = parser.parse_args()
# genetic_file = args.genetic_file
popgwas = args.popgwas
sibgwasperm = args.sibgwasperm
outdir = args.outdir
outlabel = args.outlabel
ancdata = args.ancdata
chrom = args.chrom
pos = args.pos
snpid = args.snpid
pval = args.P

alt_allele = args.alt_allele
snpset = args.snpset
sib_beta = args.sib_beta
standard_beta = args.standard_beta
log10p = args.log10p
logp_col = args.logp_col

class make_input_files(object):

	def __init__(self, outdir, outlabel, pop_gwas_file, sib_gwas_perm_file, anc_data, 
		chr_label, pos_label, snpid, alt_allele, standard_beta, sib_beta, log10p, logp_col, snpset = None): # genetic_file, 

		# Check for all files before beginning the analysis
		file_list = [pop_gwas_file, sib_gwas_perm_file, anc_data] #, genetic_file]
		for f in file_list:
			if os.path.exists(f) == False:
				raise FileNotFoundError(f'{f} does not exist or could not be opened.')

		self.pop_gwas_file_name = pop_gwas_file
		self.sib_gwas_file_name = sib_gwas_perm_file
		self.anc_data_file_name = anc_data
		self.anc_data = pd.read_csv(anc_data, sep = '\t', compression='infer')
		# self.genetic_file = genetic_file
		self.standard_beta = standard_beta
		self.sib_beta = sib_beta
		self.log10p = log10p
		self.outdir = outdir
		self.outlabel = outlabel
		self.chr_label = chr_label
		self.pos_label = pos_label
		self.snpid = snpid
		self.alt_allele = alt_allele
		self.logp_col = logp_col

		print(f"[args] Population GWAS file: {self.pop_gwas_file_name}")
		print(f"[args] Sibling GWAS file: {self.sib_gwas_file_name}")
		print(f"[args] Reference allele file: {self.anc_data_file_name}")
		# print(f"[args] Genotype file: {self.genetic_file}")
		print(f"[args] Output directory: {self.outdir}")
		print(f"[args] Output file stem: {self.outlabel}")
		print(f"[args] Population GWAS effect column: {self.standard_beta}")
		print(f"[args] Sibling GWAS effect column: {self.sib_beta}")
		if self.log10p == True:
			print("[args] --p-is-log: ON")
			print(f"[args] log10p column: {self.logp_col}")

		print(f"\nMunge sumstats began at {datetime.now()}")

		if '.gz' in self.pop_gwas_file_name:
			self.pop_gwas = pd.read_csv(self.pop_gwas_file_name, sep='\t', compression = 'gzip')
		else:
			self.pop_gwas = pd.read_csv(self.pop_gwas_file_name, sep='\t')
		# print(self.pop_gwas.head(2))
		# exit()
		
		if '.gz' in self.sib_gwas_file_name:
			self.sib_gwas = pd.read_csv(self.sib_gwas_file_name, delim_whitespace = True, compression = 'gzip')
		else:
			self.sib_gwas = pd.read_csv(self.sib_gwas_file_name, delim_whitespace = True)
		
		self.pop_gwas = self.pop_gwas.rename(columns = {snpid:'SNP'})
		pop_gwas_total = len(self.pop_gwas)
		print(f"{pop_gwas_total} variants loaded from population GWAS.")

		self.sib_gwas = self.sib_gwas.rename(columns = {snpid:'SNP'})
		sib_gwas_total = len(self.sib_gwas)
		print(f"{sib_gwas_total} variants loaded from sibling GWAS.")

		self.sib_gwas = self.sib_gwas.drop_duplicates(subset = ['SNP'])
		# print(f"LINE 83 sibgwas: {self.sib_gwas.shape}")
		self.pop_gwas = self.pop_gwas.drop_duplicates(subset = ['SNP'])
		# print(f"LINE 85 popgwas: {self.pop_gwas.shape}")

		self.sib_gwas = self.sib_gwas.dropna()
		self.pop_gwas = self.pop_gwas.dropna()
		print(f"{len(self.pop_gwas) - pop_gwas_total} variants removed as NA or duplicates from population GWAS.")
		print(f"{len(self.sib_gwas) - sib_gwas_total} variants removed as NA or duplicates from sibling GWAS.")

		if snpset:
			snpset = pd.read_csv(snpset, header = None)
			snpset.columns = ['SNP']
			print(f"{len(snpset)} variants loaded from preselected SNP file.")
		else:
			snpset = pd.DataFrame(np.zeros((1,2)))
		
		self.check_shared_snps(snpset)

	# def extract_and_clump(self):
		
	# 	if os.path.isfile(self.outdir + '/' + self.outlabel + '.PRS.clumps.snp'):
	# 		pass
	# 	else:
	# 		os.system('PRSice_linux --base ' + outdir + '/' + self.outlabel + '.support.overlap.linear --a1 ' + self.alt_allele + ' --beta --pvalue ' + pval + ' --snp SNP --target ' + self.genetic_file + ' --bar-levels 1 --fastscore --no-regress --clump-r2 0.1 --clump-kb 100 --print-snp --out ' + self.outdir + '/' + self.outlabel + '.PRS.clumps')
	# 		os.system("""awk 'NR!=1 {print $2}' """ + self.outdir + """/""" + self.outlabtl + """.PRS.clumps.snp > """ + self.outdir +"""/""" + self.outlabel + """.clumped.snps.txt""")

	def check_alt_consensus(self):
		# self.anc_data['SNP'] = self.anc_data['CHR'].astype(str) + ':' + self.anc_data['POS'].astype(str)
		relevant_anc = self.anc_data[['SNP','alt.allele']]
		print(f"{len(relevant_anc)} variants loaded from ancestral data file.")

		#check effects for to make sure that the alternative allele matches between the summary statistics and the target cohort
		self.pop_gwas = self.pop_gwas.merge(relevant_anc, on = 'SNP', how = 'inner')
		self.sib_gwas = self.sib_gwas.merge(relevant_anc, on = 'SNP', how = 'inner')
		# print(f"LINE 111 pop x relevant_anc: {self.pop_gwas.shape}")
		matcher = self.sib_gwas[['SNP']].merge(self.pop_gwas[['SNP']], on = 'SNP', how = 'inner')
		# print(f"LINE 113 sib x pop x relevant_anc: {matcher.shape}")
		self.pop_gwas = self.pop_gwas[self.pop_gwas['SNP'].isin(matcher['SNP'])].drop_duplicates(subset = ['SNP'])
		# print(f"LINE 115 pop x matcher: {self.pop_gwas.shape}")
		self.sib_gwas = self.sib_gwas[self.sib_gwas['SNP'].isin(matcher['SNP'])].drop_duplicates(subset = ['SNP'])
		# print(f"LINE 117 sib x matcher: {self.sib_gwas.shape}")
		self.pop_gwas = self.pop_gwas.reset_index(drop=True)
		self.sib_gwas = self.sib_gwas.reset_index(drop=True)

		self.pop_gwas = self.pop_gwas[self.pop_gwas['SNP'].isin(matcher['SNP'])].drop_duplicates(subset = ['SNP'])
		self.sib_gwas = self.sib_gwas[self.sib_gwas['SNP'].isin(matcher['SNP'])].drop_duplicates(subset = ['SNP'])
		# print(f"LINE 123 pop nodup: {self.pop_gwas.shape}")
		# print(f"LINE 124 sib nodup: {self.sib_gwas.shape}")

		self.pop_gwas['effect_matches'] = (self.pop_gwas[self.alt_allele]==self.pop_gwas['alt.allele']).astype(int)
		self.pop_gwas['effect_matches'] = 2*(self.pop_gwas['effect_matches']) - 1
		self.pop_gwas['effect_matches'] = self.pop_gwas['effect_matches'].astype(float)
		self.pop_gwas['beta.altconsensus'] = self.pop_gwas[self.standard_beta].astype(float) * self.pop_gwas['effect_matches']
		self.sib_gwas['effect_matches'] = (self.sib_gwas[self.alt_allele]==self.sib_gwas['alt.allele']).astype(int)
		self.sib_gwas['effect_matches'] = 2*(self.sib_gwas['effect_matches']) - 1
		self.sib_gwas['effect_matches'] = self.sib_gwas['effect_matches']
		self.sib_gwas['beta.altconsensus'] = self.sib_gwas[self.sib_beta].astype(float) * self.sib_gwas['effect_matches']

	def check_shared_snps(self, preselected_snp_ids):
		self.check_alt_consensus()
		if os.path.isfile(self.outdir + '/' + self.outlabel + '.support.overlap.linear') and os.path.isfile(self.outdir + '/' + self.outlabel + '.clumped.snps.txt') and preselected_snp_ids.shape[0] == 1:

			clumped_snps = pd.read_csv(self.outdir + '/clumped.snps.txt',sep = '\t', header = None)
			clumped_snps.columns = ['SNP']
			
			self.pop_gwas['ID'] = self.pop_gwas[self.chr_label].astype(str) + ':' + self.pop_gwas['SNP'].astype(str)
			self.pop_gwas.to_csv(self.outdir + '/' + self.outlabel + '.support.overlap.linear', index = False, sep = '\t')
			self.pop_gwas = self.pop_gwas[self.pop_gwas['SNP'].isin(clumped_snps['SNP'].tolist())].drop_duplicates()
			self.sib_gwas['ID'] = self.sib_gwas['CHR'].astype(str) + ':' + self.sib_gwas['SNP'].astype(str)
			# self.sib_gwas['SNP'] = self.sib_gwas['CHR'].astype(str) + ':' + self.sib_gwas['BP'].astype(str)
			self.sib_gwas = self.sib_gwas[self.sib_gwas['SNP'].isin(self.pop_gwas['SNP'].tolist())].drop_duplicates()
			self.pop_gwas = self.pop_gwas.reset_index(drop=True)
			self.sib_gwas = self.sib_gwas.reset_index(drop=True)
			# print(f"LINE 148 pop nodup: {self.pop_gwas.shape}")
			# print(f"LINE 149 sib nodup: {self.sib_gwas.shape}")
			# sys.exit()

		elif preselected_snp_ids.shape[0] > 1:
			# print(f"LINE 153 preselected_snp_ids: {preselected_snp_ids.shape}")
			val = preselected_snp_ids.drop_duplicates(subset = ["SNP"]).shape
			# print(f"LINE 154 preselected_snp_ids: {val}")
			# print(self.pop_gwas.filter(regex="1:717587")) ############### ADDED TO SOLVE BUG WITH MERGE ###############
			# print(self.pop_gwas.filter(regex="1:1194638")) ############### ADDED TO SOLVE BUG WITH MERGE ###############
			# preselected_snp_ids['SNP'] = preselected_snp_ids['SNP'].astype(str) ############### ADDED TO SOLVE BUG WITH MERGE ###############
			# self.pop_gwas['SNP'] = self.pop_gwas['SNP'].astype(str) ############### ADDED TO SOLVE BUG WITH MERGE ###############
			self.pop_gwas = self.pop_gwas.drop_duplicates(subset = ['SNP'])
			# print(f"LINE 164 pop nodup: {self.pop_gwas.shape}")
			self.pop_gwas = self.pop_gwas.merge(preselected_snp_ids, on = 'SNP', how = 'inner')
			# print(f"LINE 165 pop x pssnp: {self.pop_gwas.shape}")
			# missing = preselected_snp_ids[~preselected_snp_ids["SNP"].isin(self.pop_gwas["SNP"])]
			# print(missing.head(10))


			# self.pop_gwas = pd.merge(self.pop_gwas, preselected_snp_ids, on = 'SNP', how = 'inner') ############### ADDED TO SOLVE BUG WITH MERGE ###############
			self.pop_gwas.to_csv(self.outdir + '/' + self.outlabel + '.support.overlap.linear', index = False, sep = '\t')
			self.sib_gwas = self.sib_gwas.merge(preselected_snp_ids, on = 'SNP', how = 'inner')
			# print(f"LINE 159 sib x pssnp: {self.sib_gwas.shape}")
			merged = self.pop_gwas.merge(self.sib_gwas, on = 'SNP', how = 'inner')
			# print(f"LINE 161 pop x snp x pssnp: {merged.shape}")
			self.pop_gwas = self.pop_gwas.loc[self.pop_gwas['SNP'].isin(merged['SNP'])]
			self.sib_gwas = self.sib_gwas[self.sib_gwas['SNP'].isin(merged['SNP'])]
			self.pop_gwas = self.pop_gwas.reset_index(drop=True)
			self.sib_gwas = self.sib_gwas.reset_index(drop=True)

			self.clump_gwas = self.pop_gwas[[self.chr_label, 'SNP', self.pos_label, alt_allele, self.standard_beta, pval]]
			self.clump_gwas = self.clump_gwas.drop_duplicates(subset = 'SNP')
			self.clump_gwas.to_csv(self.outdir + '/' + self.outlabel + '.support.overlap.linear', index = False, sep = '\t')
			# print(f"LINE 170 pop nodup: {self.pop_gwas.shape}")
			# print(f"LINE 171 sib nodup: {self.sib_gwas.shape}")
			# sys.exit()


		else:
			#first make sure that the alternative allele is set to be the same as in the 1kg data
			#now that we have fixed that let's merge them together, using the chromosome rsid combo label to
			#account for any potential repetitions of rsIDs on different chromosomes
			#take the shared SNPs, relabel them, and then write out the list
			shared_snps = self.pop_gwas[['SNP']].merge(self.sib_gwas[['SNP']], on = 'SNP', how = 'inner')
			consensus = shared_snps[['SNP']].merge(self.anc_data[['SNP']], on = 'SNP', how = 'inner')

			#now that we have the overlap between snps sets take the corresponding summary statistics from 
			#the standard GWAS and clump them agnostic of p-value
			self.pop_gwas = self.pop_gwas.set_index('SNP').loc[consensus['SNP'].tolist()]
			self.pop_gwas = self.pop_gwas.rename(columns ={'TEST.1':'STAT','OBS_CT':'NIND'})
			self.pop_gwas = self.pop_gwas.reset_index()
			self.clump_gwas = self.pop_gwas[[self.chr_label,'SNP','BP', 'A1', 'BETA',pval]]
			self.clump_gwas = self.clump_gwas.drop_duplicates(subset = 'SNP')
			self.clump_gwas.to_csv(self.outdir + '/' + self.outlabel + '.support.overlap.linear', index = False, sep = '\t')
			#extract the consensus from the provided 1kg file and clump them agnostic to p-value
			# self.extract_and_clump()

			#read in the resulting SNPs from each clump
			clumped_snps = pd.read_csv(self.outdir + '/' + self.outlabel + '.clumped.snps.txt',sep = '\t', header = None)
			clumped_snps.columns = ['SNP']
			self.pop_gwas = self.pop_gwas.merge(clumped_snps, on = 'SNP', how = 'inner').drop_duplicates()
			self.sib_gwas = self.sib_gwas.reset_index(drop = True)
			self.sib_gwas = self.sib_gwas[self.sib_gwas['SNP'].isin(self.pop_gwas['SNP'].tolist())]

					
		self.pop_gwas = self.pop_gwas.rename(columns={'SE':'se'})
		
		if '#CHROM' in self.pop_gwas.columns:
			self.pop_gwas = self.pop_gwas.rename(columns={'#CHROM':'CHR'})
		
		if self.log10p:
			self.pop_gwas['P'] = 10**(-1.0*self.pop_gwas[self.logp_col])

		if len(self.pop_gwas) < 1:
			raise RuntimeError("ERROR: No SNPs passed filtering and munging for population GWAS.")
		else:
			print(f"{len(self.pop_gwas)} variants remain in population GWAS after filtering and merging.")
		if len(self.sib_gwas) < 1:
			raise RuntimeError("ERROR: No SNPs passed filtering and munging for sibling GWAS.")
		else:
			print(f"{len(self.sib_gwas)} variants remain in sibling GWAS after filtering and merging.")

		self.pop_gwas.to_csv(self.outdir + '/' + self.outlabel + '.standard.preproc.txt', sep = '\t', index = False)
		siblabel = self.sib_gwas_file_name.replace('.gz','')
		siblabel = self.sib_gwas_file_name.replace('.txt.gz','')
			
		self.sib_gwas.to_csv(self.outdir + '/' + self.outlabel + '.sib.preproc.txt', sep = '\t', index = False)

if __name__ == '__main__':
	x = make_input_files(outdir, outlabel, popgwas, sibgwasperm, ancdata, chrom, pos, snpid,
	 alt_allele, standard_beta, sib_beta, log10p, logp_col, snpset) # genetic_file, 
	print(f'Munge complete at {datetime.now()}.\n')