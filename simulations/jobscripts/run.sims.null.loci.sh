#!/bin/bash
#SBATCH -J test.sim
#SBATCH -o slurm_outputs/test.sim.o%j
#SBATCH -e test.sim.o%j
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH --mail-user=samuel.smith@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH -A OTH21148
#SBATCH --array=1-5
#These simulations represent the null architecture with varying direct and and environmental variance
#maf linked genetic architecture
source activate sad3.9
# echo $SLURM_ARRAY_TASK_ID

#env-var equal 0.99
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
	# python DataDriver.py --neutral-architecture --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1.0.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.0.001.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-05.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-08.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	
	# python DataDriver.py --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1.0.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.0.001.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-05.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.01 --direct-variance 0.99 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-08.npcs.20.env.var.0.01.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	
fi
#env-var equal 0.9
if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
	# python DataDriver.py --neutral-architecture --nsim 63 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1.0.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.0.001.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-05.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-08.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	
	# python DataDriver.py --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1.0.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.0.001.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-05.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.1 --direct-variance 0.9 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-08.npcs.20.env.var.0.1.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	
fi

#env-var equal 0.8
if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then	
	# python DataDriver.py --neutral-architecture --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1.0.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.0.001.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-05.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-08.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

	# python DataDriver.py --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1.0.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.0.001.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-05.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.2 --direct-variance 0.8 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-08.npcs.20.env.var.0.2.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

fi

if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then
	# python DataDriver.py --neutral-architecture --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1.0.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.0.001.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-05.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-08.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

	# python DataDriver.py --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1.0.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.0.001.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-05.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.5 --direct-variance 0.5 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-08.npcs.20.env.var.0.5.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

fi

if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then
	# python DataDriver.py --neutral-architecture --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1.0.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.0.001.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-05.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	python DataDriver.py --neutral-architecture --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.neutral.architecture.pval.1e-08.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

	# python DataDriver.py --nsim 100 --ascp 1.0 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1.0.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.0.001.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-05.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400
	# python DataDriver.py --nsim 100 --ascp 0.00000001 --ngwaspcs 20 --n-trait-loci 400 --env-variance 0.8 --direct-variance 0.2 --indirect-variance 0 --direct-indirect-covar 0 --env-covariance 0.5 --n-sib-gwas 20000 --n-standard-gwas 50000 --pc-lower-bound 100 --out-label sib.nosad.maf.architecture.pval.1e-08.npcs.20.env.var.0.8.env.cov.0.5.nsib.20000.nstandard.50000.ntraitloci.400

fi
