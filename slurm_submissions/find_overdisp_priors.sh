#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=4G 
#SBATCH -o update_estimates_2-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=0-56

module purge
module load R
cd //dfs6/pub/igoldst1/CA_rt_estimates


Rscript scripts/find_overdisp_priors.R $SLURM_ARRAY_TASK_ID
