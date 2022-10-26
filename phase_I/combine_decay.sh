#!/bin/bash


#SBATCH -A ALLISONS_LAB             ## account to charge
#SBATCH -p standard                 ## run on the standard partition
#SBATCH -N 1                        ## run on a single node
#SBATCH --mem 3GB                   ## memory requested
#SBATCH --ntasks-per-node=1         ## number of tasks to launch per node
#SBATCH --error=slurm-%J.err        ## write errors in slurm-<jobID>.err file
#SBATCH --mail-type=end             ## send email when the job ends
#SBATCH --mail-user=eabs@uci.edu    ## use this email address

module purge
hostname
module load python/3.10.2
cd src
python combine_decay.py
