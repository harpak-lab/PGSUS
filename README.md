# PGSUS software for the estimation of SAD variance in a polygenic score

The Partitioning Genetic Scores Using Siblings (PGSUS, pron. "Pegasus") method was developed to decompose the variance of a polygenic score (PGS) into its different sources of variation. Namely, the PGSUS method uses paried population level GWAS statistics and sibling GWAS statistics to identify the proportion of variance in a PGS due to direct effects and variance due to Stratification, Assortative mating, and indirect parental (Dynastic) effects or ``SAD effects''. The PGSUS method is written in Python and run in conda environment described and discussed below. Notably, application of PGSUS helps in the identification and interpretation of the performance of a PGS with respect to a particlar target cohort. It achieves this through a two-step decomposition of variance in a PGS that is due to SAD effects that is shared across all of the PCs of the genotype matrix of the prediction sample and the SAD variance that is correlated with a particular axis of stratification, as captured by the prediction sample PCs.

## Directory Contents
In addition to the files below, files containing the summary statistics, preprocessing steps, and actual commands for each analysis can be found on the [Harpak Lab Website Data Tab](https://www.harpaklab.com/data). 

`simulations` contains the scripts used to generate simulated summary statistics and perform the PGSUS decomposition while varying parameters of interst. Further discussion of the implementation and different flags can be found in this directory. 

## Formatting of summary statistics
The first step in applying PGSUS is formatting the standard and sibling GWAS summary statistics properly. This is achieved using the `munge_sumstats.py` script and the follwoing set of flags. In addition, **please be sure to download the necessary support files** from the link above. An example command, and enumeration of each flag can be seen below. 

```python 
python munge_sumstats.py --pop-gwas-file standard_gwas_height.linear.gz
--sib-perm-file sibling_gwas_height.txt.gz
--preselected-snps index_snps.txt
--outdir myresults/
--outlabel pgsus_height_decomposition
--chr CHR
--snp-id ID
 ```

`--pop-gwas-file` 


