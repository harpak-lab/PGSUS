# PGSUS software for the estimation of SAD variance in a polygenic score

The Partitioning Genetic Scores Using Siblings (PGSUS, pron. "Pegasus") method was developed to decompose the variance of a polygenic score (PGS) into its different sources of variation. Namely, the PGSUS method uses paired population level GWAS statistics and sibling GWAS statistics to identify the proportion of variance in a PGS due to direct effects and variance due to Stratification, Assortative mating, and indirect parental (Dynastic) effects or "SAD effects". The PGSUS method is written in Python and run in conda environment described and discussed below. Notably, application of PGSUS helps in the identification and interpretation of the performance of a PGS with respect to a particlar target cohort. It achieves this through a two-step decomposition of variance in a PGS that is due to SAD effects that is shared across all of the PCs of the genotype matrix of the prediction sample and the SAD variance that is correlated with a particular axis of stratification, as captured by the prediction sample PCs.

## Directory Contents
In addition to the files below, files containing the summary statistics, preprocessing steps, and actual commands for each analysis can be found on the [Harpak Lab Website Data Tab](https://www.harpaklab.com/data). 

`simulations` contains the scripts used to generate simulated summary statistics and perform the PGSUS decomposition while varying parameters of interst. Further discussion of the implementation and different flags can be found in this directory. 

## Formatting of summary statistics
The first step in applying PGSUS is formatting the standard and sibling GWAS summary statistics properly. This is achieved using the `munge_sumstats.py` script and the follwoing set of flags. In addition, **please be sure to download the necessary support files** from the link above. An example command where the index variants are already known and enumeration of each flag can be seen below. 

```python 
python munge_sumstats.py --pop-gwas-file standard_gwas_height.linear.gz
--sib-perm-file sibling_gwas_height.txt.gz
--preselected-snps index_snps.txt
--outdir myresults/
--outlabel pgsus_height_decomposition
--snp-id ID
 ```

The possible flags that can be used to munge different input file formats are enumerated below:
`--pop-gwas-file` a file containing summary statistics from a population GWAS. In this example, a compressed file from [plink](https://www.cog-genomics.org/plink/1.9/assoc#linear)'s implementation is used. 

`--sib-perm-file` a file containing summary statistics from a sibling GWAS. 

`--preselected-snps` a file with a single column containing the set of SNP IDs to be used in the PGSUS decomposition.

`--outdir` directory path where results should be written. 

`--outlabel` file prefix for output files produced during the data munging. 

`--snp-id` label of then column containing SNP IDs formatted as chromsome:basepair. Defaults to "SNP". 

`--chr` label of the column containing the chromosome of each SNP. Default is "CHR". 

`--bfile` path to the plink formatted bed file with the target sample. Only required if the clumping of the loci has not been previously been performed and identified with the preselected SNPs flag.

`--anc-data` file containing ancestral and derived alleles for each locus. Defaults to the file created using the 1000 Genomes data contained in the support files directory. 

`--pos` label of the base pair position of each locus. Defaults to "POS".

`--standard-beta` label of the column in the population GWAS file that contains the effect size estimates. Defaults to "BETA".

`--sib-beta` label of the column in the sibling GWAS file that contains the effect size estimates. Defaults to "BETA".

`--pval` label of the column containing the population GWAS association p-values. Defaults to "P".

`--a1` label of the column containing the alternate allele counted in the GWAS for each locus. Defaults to "A1". 

`--p-is-log` presence of this flag indicates that p-values have been negative log 10 transformed. 

`--logp-col` label of the column containing the negative log 10 transformed p-values. Only necessary if the p-value has undergone transformation. 

Once the script is run, there should be two files produced both with specified outlabel above: one with the suffix `sib.preproc.txt` and another with the the suffix `standard.preproc.txt` each containing the formatted summary statistics from the population and sibling GWAS. Note that there will be a new column entitled `beta.altconsensus` that contains the same effect size estimates, but some will have been multiplied by -1 in order to ensure the correct polarity of the effect allele. 


## Estimation of SAD variance using the PGSUS software

With the preprocessed data in hand the PGSUS software can now be run. There are a number of possible flags intended to help incorporate different sources of data. An example command is given below and each possible flag is described as well. 

```python
python pgsus.py --genetic-file 1KG.data.bed
--pop-gwas myresults/pegasus_height_decomposition.standard.preproc.txt
--sib-gwas myresults/pegasus_height_decomposition.sib.preproc.txt
--chrom-pos SNP
--pvalue 1
--pval-col P
--pop-effect beta.altconsensus
--pop-se se
--sib-effect beta.altconsensus
--sib-se EMP_SE
--ascertainment-set gwas
--outfile-label plink.wc.1kg.all.sps23.height.aperm.1K.to.1M
--out sps_23_stats_v2/height/
--chrom CHR
--pos BP
--permutation-test
```

`--genetic-file` the plink formatted gentoype file for the prediction sample. This should contain the same SNPs that are in the ".preproc.txt" files. If not, the overlap will be found. 

`--pop-gwas-file` the preprocessed file containing the population GWAS summary statistics. 

`--sib-gwas-file` the preprocessed file containing the sibling GWAS summary statistics. 

`--pvalue` the label of the column containing the population GWAS association pvalues in the `pop-gwas-file`. 

`--ascertainment-set` either needs to be set to "gwas" or "sibs". Determines which GWAS to use for SNP selection given a particular threshold. Default is "gwas".

`--chrom` label of the column containing the chromosome of each SNP. Default is "CHR". 

`--pos` label of the base pair position of each locus. Defaults to "POS".

`--pop-effect` label of the column in the population GWAS file that contains the effect size estimates. Defaults to "BETA".

`--pop-se` label of the column in the population GWAS file that contains the standard errors. Defaults to "se". 

`--sib-effect` label of the column in the sibling GWAS file that contains the effect size estimates. Defaults to "BETA".

`--sib-se` label of the column in the sibling GWAS file that contains the standard errors. Defaults to "se". 

`--pval-col` label of the column containing the population GWAS association p-values. Defaults to "P".

`--nboots` the number of bootstraps to perform in estimating the standard errror of the isotropic inflation factor using the Deming regression framework.

`--eigenvals` file that contains the eigenvalues of the prediction sample genotype matrix saved as an ".npy" formatted file. If this flag is absent then a numpy implementation of PCA will be berformed on the plink binary file specified using `bfile`. 

`--eigenvecs` file that contains the eigenvectors of the prediction sample genotype matrix saved as an ".npy" formatted file.

`--permutation-test` the presence of this flag indicates that the PC-wise permutation test should be performed. If absent, only the isotropic inflation factor will be estimated. 

`--perm-pcs` the number of top PCs to be tested for the PC-wise permutation procedure. Defaults to the first 100 PCs. 

`--nperm` the number of permutations to perform for each PC-wise decomposition in constructing the permutation based null. Defaults to 1000 permutations. 

`--outfile-label` file prefix for output files produced during the data munging. 

`--out` directory path where results should be written. 

