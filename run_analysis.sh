#!/bin/bash

#SBATCH --job-name=run_traj_analysis
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-type=ALL

### pre-configuration ###

data_prefix=Kriegstein
run_path=$(pwd)
data_path=/gpfs/gibbs/pi/gerstein/pse5/Hybrid_cell_scheme/
use_coding=FALSE

mkdir -p ${run_path}/data
mkdir -p ${run_path}/data/gam
mkdir -p ${run_path}/data/gam/${data_prefix}_gam
mkdir -p ${run_path}/results
mkdir -p ${run_path}/ds_analysis


### load conda environment ###

module load miniconda
conda init bash
source activate /gpfs/ysm/project/gerstein/ah2426/conda_envs/R4.2.0


### run pipeline

# step 1: make seurat object & run trajectory analysis
Rscript ${run_path}/1_makeSeurat_trajAnalysis.R ${data_prefix} ${run_path} ${data_path} ${use_coding}

# step 2: make the design matrix for fitGAM
Rscript ${run_path}/2_makeGAM.R ${data_prefix} ${run_path}

# step 3: fit GAM
Rscript ${run_path}/3_fitGAM.R ${data_prefix} ${run_path}

# NOTE: if want to fitGAM in parallel, need to do it outside of this file 
