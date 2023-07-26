#!/bin/bash
#SBATCH -J test.sim
#SBATCH -o test.sim.o%j
#SBATCH -e test.sim.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH --mail-user=samuel.smith@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --array=1-20


pval_array=(0.00000001  0.0000001  0.000001 0.00001 0.0001 0.001 0.01 0.1 1)
lb_array=(10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200)

# for p_cutoff in "${pval_array[@]}"
# do
# 	python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp $p_cutoff --pc-lower-bound 100
# done


for p_cutoff in "${pval_array[@]}"
do
	python DataDriver.py --nsim 100 --ngwaspcs 5 --ascp $p_cutoff --pc-lower-bound ${lb_array[$SLURM_ARRAY_TASK_ID]}
done

