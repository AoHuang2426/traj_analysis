#!/bin/bash

#SBATCH --job-name=after_fitGAM_parallel
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=ALL

#########################
### pre-configuration ###
#########################

data_prefix=Ma_Sestan
run_path=$(pwd)

#########################


### load conda environment ###

module load miniconda
conda init bash
source activate /gpfs/ysm/project/gerstein/ah2426/conda_envs/R4.2.0


### the last step of the pipeline

Rscript ${run_path}/parallel_run/2.2_combineGAM.R ${data_prefix} ${run_path}
