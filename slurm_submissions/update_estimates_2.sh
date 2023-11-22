#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 06:00:00   ## 4 hr run time limit
#SBATCH --mem=4G 
#SBATCH -o update_estimates_2-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=0-50,52-56
module purge
module load R
cd //pub/igoldst1/CA_rt_estimates

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/update_estimates_3.sh
fi

Rscript scripts/fit_estimgamma.R $SLURM_ARRAY_TASK_ID
