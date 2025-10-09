import pandas as pd 
import numpy as np 
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# for trait in ['pack_years_smoking','height']:
# 	standard = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.standard.preproc.txt', sep = '\t')
# 	sib = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.sib.preproc.txt', sep = '\t')
# 	sib['T_STAT'] = sib['BETA']/sib['EMP_SE']

# 	single_snp_clumps = pd.read_csv('../../6.empirical_estimation/omit_singleton_clumps/' + trait + '.clumped.nosingleton.ids.txt',sep = '\t',header = None)
# 	single_snp_clumps.columns = ['SNP']
# 	colormap = {True:'#377eb8', False:'#ff7f00'}
# 	fig, ax = plt.subplots(1,4,figsize = (14,4))

# 	sib_prop = sib[(sib['T_STAT'] > -1) & (sib['T_STAT'] < 1)].shape[0]/float(sib.shape[0])
# 	standard_prop = standard[(standard['T_STAT'] > -1) & (standard['T_STAT'] < 1)].shape[0]/float(standard.shape[0])
# 	ax[0].hist(standard['T_STAT'])
# 	ax[0].title.set_text('Standard Z scores\n' + str(round(standard_prop,4)*100) + r'% $-1<Z<1$')
# 	ax[0].set_xlabel('Z scores')
# 	ax[0].set_ylabel('Number of SNPs')

# 	ax[1].hist(sib['T_STAT'])
# 	ax[1].title.set_text('Sibling Z scores\n' + str(round(sib_prop,4)*100) + r'% $-1<Z<1$')
# 	ax[1].set_xlabel('Z scores')
# 	ax[1].set_ylabel('Number of SNPs')

# 	sib['standard_T_STAT'] = standard['T_STAT']
# 	ax[2].scatter(sib['T_STAT'],sib['standard_T_STAT'], alpha = 0.5)
# 	ax[2].title.set_text('Sibling v. Standard\nZ scores')
# 	ax[2].set_xlabel('Sibling Z')
# 	ax[2].set_ylabel('Standard Z')


# 	sib['single_clump'] = sib['SNP'].isin(single_snp_clumps['SNP'])
# 	sib['single_clump'] = sib['single_clump'].map(colormap)
# 	multi_clumps = sib[sib['single_clump'] == '#377eb8']
# 	ax[3].scatter(multi_clumps['T_STAT'],multi_clumps['standard_T_STAT'], alpha = 0.1)
# 	single_clumps = sib[sib['single_clump'] == '#ff7f00']

# 	single_clumps_prop = single_clumps[(single_clumps['T_STAT'] > -1) & (single_clumps['T_STAT'] < 1)].shape[0]/float(single_clumps.shape[0])

# 	ax[3].scatter(single_clumps['T_STAT'],single_clumps['standard_T_STAT'], alpha = 0.9, color = 'orange')
# 	ax[3].title.set_text('Sibling v. Standard\nZ scores single clumps\n' + str(round(single_clumps_prop,4)*100) + r'% $-1<Z<1$')
# 	ax[3].set_xlabel('Sibling Z')
# 	ax[3].set_ylabel('Standard Z')

# 	sns.despine()
# 	plt.tight_layout()
# 	plt.savefig('se.check.' + trait + '.png')
# 	plt.clf()

sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
		'hand_grip_strength','height','hip_circ','household_income','neuroticism_score',
		'overall_health','pack_years_smoking','pulse_rate','skin_color','waist_circ','years_schooling']
thresholds = [1,1e-3,1e-5,1e-8]

# for trait in sps_traits:
# 	fig, ax = plt.subplots(1,4,figsize = (14,4))
# 	standard = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.standard.preproc.txt', sep = '\t')
# 	sib = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.sib.preproc.txt', sep = '\t')
# 	standard['sib_T_STAT'] = sib['BETA']/sib['EMP_SE']

# 	for j,k in enumerate(thresholds):
# 		temp = standard[standard['P'] < k].dropna()
# 		print(trait,np.min(temp['sib_T_STAT']), np.max(temp['sib_T_STAT']))
# 		ax[j].hist(temp['sib_T_STAT'])
# 		sib_prop = temp[(temp['sib_T_STAT'] > -1) & (temp['sib_T_STAT'] < 1)].shape[0]/float(temp.shape[0])
# 		ax[j].title.set_text('Sibling Z scores\nStandard ascertainment\n' + str(round(sib_prop,4)*100) + r'% $-1<Z<1$, nSNP = ' + str(temp.shape[0]) )

# 	sns.despine()
# 	plt.tight_layout()
# 	plt.savefig('sib.z.threhsolds.' + trait + '.pdf')


# outdf = pd.DataFrame(np.zeros((len(sps_traits),thresholds)), index = sps_traits, columns = ['p_passed'])

# for trait in sps_traits:
# 	sib = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.sib.preproc.txt', sep = '\t')
# 	standard = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.wc.' + trait + '.aperm.1K.to.1M.standard.preproc.txt', sep = '\t')
# 	outdf.loc[trait,'p_passed'] = sib[sib['EMP1'] < 0.000005].shape[0]

# outdf.to_csv('nSNP.aperm.max.txt',sep = '\t')

# napermmax = pd.read_csv('nSNP.aperm.max.txt', sep = '\t')
# alpha_snps = pd.read_csv('../cache/alpha_matrices/plink.wc.1kg.all.sps23.aperm.1K.to.1M.v2.nsnp.mat.txt', sep = '\t')

# fig, ax = plt.subplots(1,1,figsize=(4,4))
# print(napermmax)

# ax.scatter(alpha_snps['1e-08'],napermmax['p_passed'])
# ax.set_xlabel(r'nSNPs in PGS at $p<10^{-8}$')
# ax.set_ylabel('nSNPs with permutation\n' + r'$p<5\times10^{-6}$')
# print(np.corrcoef(alpha_snps['1e-08'],napermmax['p_passed'])**2)
# plt.tight_layout()
# sns.despine()
# plt.savefig('perm.v.alpha.1e-8.snps.pdf')
# plt.clf()
import scipy.stats as stats

sps_traits = ['alcohol_intake_freq','birth_weight','bmi','dbp','fvc','hair_color',
		'hand_grip_strength','height','hip_circ','household_income',
		'pulse_rate','skin_color','waist_circ','years_schooling']

sps_traits = ['height']
for trait in sps_traits:
	ah = pd.read_csv('preproc_mh20_height_sib.tsv',sep = '\t')
	ah = ah.rename(columns={'chrom.pos':'SNP'})
	ah['norm_se'] = np.abs(ah['BETA']/stats.norm.ppf(ah['RAW_P'].astype(float)/2))
	# sys.exit()
	# sib_normal = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.qfam.normal.theory.' + trait + '.aperm.1K.to.1M.se.v2.txt.gz', sep = '\t', compression = 'gzip')
	# sib_normal['SNP'] = sib_normal['CHR'].astype(str) + ':' + sib_normal['BP'].astype(str)
	sib_perm = pd.read_csv('../../6.empirical_estimation/sps_23_stats_v2/' + trait + '/plink.qfam.perm.' + trait + '.aperm.1K.to.1M.se.v2.txt.gz', sep = '\t', compression = 'gzip')
	# sib_normal['norm_se'] = np.abs(sib_normal['BETA']/stats.norm.ppf(sib_normal['RAW_P'].astype(float)/2))
	allsnps = sib_perm[sib_perm['SNP'].isin(ah['SNP'].tolist())]['SNP'].tolist()
	# sib_normal = sib_normal.set_index('SNP')
	ah = ah.set_index('SNP')
	sib_perm = sib_perm.set_index('SNP')
	sib_perm = sib_perm.loc[allsnps]
	ah = ah.loc[allsnps]
	fig, ax = plt.subplots(1, 2, figsize = (8,4))
	ax[0].scatter(ah['norm_se'],sib_perm['EMP_SE'])
	ax[0].set_xlabel('Normal SE')
	ax[0].set_ylabel('Permutation SE')
	ax[1].scatter(np.negative(np.log10(ah['RAW_P'])),np.negative(np.log10(sib_perm['EMP1'])))
	ax[1].set_xlabel('Normal P')
	ax[1].set_ylabel('Permutation P')
	plt.tight_layout()
	plt.savefig(trait + '.perm.v.normal.png')
	plt.clf()




