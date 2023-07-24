#!/bin/bash
#SBATCH -J test.sim
#SBATCH -o test.sim.o%j
#SBATCH -e test.sim.o%j
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH --mail-user=samuel.smith@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --array=1-20

# python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.1 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.1 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.01 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.001 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.0001 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.00001 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.000001 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.0000001 --pc-lower-bound 100
python DataDriver.py --nsim 100 --ngwaspcs $SLURM_ARRAY_TASK_ID --ascp 0.00000001 --pc-lower-bound 100