#!/bin/bash


#SBATCH -A ALLISONS_LAB             ## account to charge
#SBATCH -N 1                        ## run on a single node
#SBATCH --mem 15GB                  ## memory requested
#SBATCH --ntasks-per-node=4         ## number of tasks to launch per node
#SBATCH --error=slurm-%J.err        ## write errors in slurm-<jobID>.err file
#SBATCH --mail-type=end             ## send email when the job ends
#SBATCH --mail-user=eabs@uci.edu    ## use this email address
#SBATCH --array=1-10

module purge
hostname
module load python/3.8.0
cd src
python dementpy.py ${SLURM_JOBID} input output ${i} ${j} ${k} 3 $SLURM_ARRAY_TASK_ID 3 
#Arg:jobID,input_file_name,output_file_name,final_litter,res,inv,enzmax,seed_index,pulse
