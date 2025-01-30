![PGSUS](figures_src/PGSUS_logo.png)
# PGSUS: Partitioning Genetic Scores Using Siblings

The Partitioning Genetic Scores Using Siblings (PGSUS, pronounced "Pegasus") method decomposes the variance of a polygenic score (PGS) into its various sources. Specifically, PGSUS uses paired population-level GWAS statistics and sibling GWAS statistics to determine the proportion of variance in a PGS attributable to direct effects and to Stratification, Assortative mating, and Indirect parental (Dynastic) effects, collectively referred to as "SAD effects." PGSUS is implemented in Python and runs in a conda environment. Applying PGSUS helps identify and interpret the performance of a PGS with respect to a particular target cohort by performing a two-step decomposition of variance in a PGS.

## Directory Contents

In addition to the files listed below, you can find files containing summary statistics, preprocessing steps, and actual commands for each analysis on the [Harpak Lab Website Data Tab](https://www.harpaklab.com/data). Additionally, the code and necessary input to generate each of the figures contained withing the manuscript can be found on the Harpak Lab Website Data Tab.

- `simulations`: contains scripts to generate simulated summary statistics and perform the PGSUS decomposition while varying parameters of interest. Detailed discussion of implementation and different flags can be found in this directory.
- `example`: contains input files needed to run the example commands in this repo.

## Creating and activating a conda environment

Before you can analyze your polygenic score, you first need to create the appropriate environment using conda. Once you have downloaded the requirements.txt file and have installed conda on your computer, the follwoing can be run to generate a conda environment with the necessary packages to apply PGSUS. 

`conda create --name pgsus --file requirements.txt`

Once the environment is successfully created activate it using

`source activate pgsus`

The PGSUS software uses pysnptools to manipulate genetic data. It is most easily installed through pip using the command below while the environment is active. 

`pip install pysnptools`


## Clumping to identify independent markers

The final prerequisite step for running PGSUS is to obtain a clumped set of SNPs that can be used to construct a PGS from independent markers. We recommend doing this using the greedy algorithm implemented by [plink1.9](https://www.cog-genomics.org/plink/1.9/postproc#clump). In our analyses we used the following command and parameter set, details for each flag can be found on the plink1.9 website. 

```plink1.9 --bfile reference_population
--clump gwas.results.linear
--clump-p1 1
--clump-p2 1
--clump-r2 0.1
--clump-kb 100
--memory 16000
--out gwas.results.clumped
```
`reference_population` is a plink binary file containing the population that will be used to estimate the LD among variants. 

`gwas.results.linear` is the output from a population GWAS containing the summary statistics that will be used to construct the PGS.

The list of index SNPs contained in the output file can then be used as the set of preselected SNPs in the following analyses. 

## Formatting Summary Statistics

The first step in applying PGSUS is properly formatting the standard and sibling GWAS summary statistics. Use the `munge_sumstats.py` script with the following flags. 

**Please download the necessary support files** from the link above before performing analyses. The support_files directory should be placed in the same directory as the pgsus.py script. 

```python 
python munge_sumstats.py --pop-gwas-file example/example.pop.stats.linear \
--anc-data example/example.alt.alleles.txt.gz \
--sib-perm-file example/example.sib.stats.linear \
--preselected-snps example/example.snp.ids.txt \
--outdir example/ \
--outlabel pgsus_height_example \
--snp-id ID 
 ```

The possible flags that can be used to munge different input file formats are enumerated below:

- `--pop-gwas-file` a file containing summary statistics from a population GWAS. In this example, a compressed file from [plink](https://www.cog-genomics.org/plink/1.9/assoc#linear)'s implementation is used. 

- `--sib-perm-file` a file containing summary statistics from a sibling GWAS. 

- `--preselected-snps` a file with a single column containing the set of SNP IDs to be used in the PGSUS decomposition.

- `--outdir` directory path where results should be written. 

- `--outlabel` file prefix for output files produced during the data munging. 

- `--snp-id` label of then column containing SNP IDs formatted as chromsome:basepair. Defaults to "SNP". 

- `--chr` label of the column containing the chromosome of each SNP. Default is "CHR". 

- `--anc-data` file containing ancestral and derived alleles for each locus. Defaults to the file created using the 1000 Genomes data contained in the support files directory. The default file can be downloaded from the dropbox and be housed in the `support_files` directory. The default will be set as `support_files/SNPalleles_1000Genomes_allsites.txt.gz`. 

- `--pos` label of the base pair position of each locus. Defaults to "POS".

- `--standard-beta` label of the column in the population GWAS file that contains the effect size estimates. Defaults to "BETA".

- `--sib-beta` label of the column in the sibling GWAS file that contains the effect size estimates. Defaults to "BETA".

- `--pval` label of the column containing the population GWAS association p-values. Defaults to "P".

- `--a1` label of the column containing the alternate allele counted in the GWAS for each locus. Defaults to "A1". 

- `--p-is-log` presence of this flag indicates that p-values have been negative log 10 transformed. 

- `--logp-col` label of the column containing the negative log 10 transformed p-values. Only necessary if the p-value has undergone transformation. 

Once the script is run, there should be two files produced both with specified outlabel above: one with the suffix `sib.preproc.txt` and another with the the suffix `standard.preproc.txt` each containing the formatted summary statistics from the population and sibling GWAS. Note that there will be a new column entitled `beta.altconsensus` that contains the same effect size estimates, but some will have been multiplied by -1 in order to ensure the correct polarity of the effect allele. 


## Estimation of SAD variance using the PGSUS software

With the preprocessed data in hand the PGSUS software can now be run. There are a number of possible flags intended to help incorporate different sources of data. An example command is given below and each possible flag is described as well. 

```python
python pgsus.py --genetic-file 1kg.example.bed \
--pop-gwas example/pgsus_height_example.standard.preproc.txt \
--sib-gwas example/pgsus_height_example.sib.preproc.txt \
--chrom-pos SNP \
--pvalue 1 \
--pval-col P \
--pop-effect beta.altconsensus \
--pop-se se \
--sib-effect beta.altconsensus \
--sib-se EMP_SE \
--ascertainment-set gwas \
--outfile-label pgsus_height_example \
--out example/  \
--chrom CHR \
--pos BP \
--permutation-test
```

Possible flags include:
- `--genetic-file` the plink formatted gentoype file for the prediction sample. This should contain the same SNPs that are in the ".preproc.txt" files. If not, the overlap will be found. 

- `--pop-gwas-file` the preprocessed file containing the population GWAS summary statistics. 

- `--sib-gwas-file` the preprocessed file containing the sibling GWAS summary statistics. 

- `--pvalue` the label of the column containing the population GWAS association pvalues in the `pop-gwas-file`. 

- `--ascertainment-set` either needs to be set to "gwas" or "sibs". Determines which GWAS to use for SNP selection given a particular threshold. Default is "gwas".

- `--chrom` label of the column containing the chromosome of each SNP. Default is "CHR". 

- `--pos` label of the base pair position of each locus. Defaults to "POS".

- `--pop-effect` label of the column in the population GWAS file that contains the effect size estimates. Defaults to "BETA".

- `--pop-se` label of the column in the population GWAS file that contains the standard errors. Defaults to "se". 

- `--sib-effect` label of the column in the sibling GWAS file that contains the effect size estimates. Defaults to "BETA".

- `--sib-se` label of the column in the sibling GWAS file that contains the standard errors. Defaults to "se". 

- `--pval-col` label of the column containing the population GWAS association p-values. Defaults to "P".

- `--nboots` the number of bootstraps to perform in estimating the standard errror of the isotropic inflation factor using the Deming regression framework.

- `--eigenvals` file that contains the eigenvalues of the prediction sample genotype matrix saved as an ".npy" formatted file. If this flag is absent then a numpy implementation of PCA will be berformed on the plink binary file specified using `bfile`. 

- `--eigenvecs` file that contains the eigenvectors of the prediction sample genotype matrix saved as an ".npy" formatted file.

- `--permutation-test` the presence of this flag indicates that the PC-wise permutation test should be performed. If absent, only the isotropic inflation factor will be estimated. 

- `--perm-pcs` the number of top PCs to be tested for the PC-wise permutation procedure. Defaults to the first 100 PCs. 

- `--nperm` the number of permutations to perform for each PC-wise decomposition in constructing the permutation based null. Defaults to 1000 permutations. 

- `--outfile-label` file prefix for output files produced during the data munging. 

- `--out` directory path where results should be written. 
