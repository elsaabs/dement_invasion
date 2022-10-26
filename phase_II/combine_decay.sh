#!/bin/bash


#SBATCH -A ALLISONS_LAB             ## account to charge
#SBATCH -N 1                        ## run on a single node
#SBATCH --ntasks-per-node=1         ## number of tasks to launch per node
#SBATCH --error=slurm-%J.err        ## write errors in slurm-<jobID>.err file
#SBATCH --mail-type=end             ## send email when the job ends
#SBATCH --mail-user=eabs@uci.edu    ## use this email address

module purge
hostname
module load python/3.8.0
cd src
python combine_decay.py
