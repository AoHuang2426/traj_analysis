#!/bin/bash

#SBATCH --job-name=run_traj_analysis
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-type=ALL

#########################
### pre-configuration ###
#########################

data_prefix=Ma_Sestan
run_path=$(pwd)
data_path=/gpfs/gibbs/pi/gerstein/pse5/Hybrid_cell_scheme/
use_coding=FALSE

run_parallel=TRUE

#########################


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
if (( ${run_parallel} == FALSE ))
then
    Rscript ${run_path}/3_fitGAM.R ${data_prefix} ${run_path}
else
    echo "Start running fitGAM in parallel"
    
    # run the scripts in ./parallel_run
    
    # step .0: make job list for fitGAM
    Rscript ${run_path}/parallel_run/2.0_makeJob.R ${data_prefix} ${run_path}
    
    # step .1: create dSQ jobs for parallel running
    module load dSQ
    dsq --job-file ${run_path}/data/gam/${data_prefix}_joblist_parallel.txt --nodes 1 --mem-per-cpu 80G --cpus-per-task 1 --time 3:00:00 --mail-type ALL
    
    mkdir -p ${run_path}/jobs-dSQ_makeGAM
    mv ${run_path}/dsq-${data_prefix}* ./jobs-dSQ_makeGAM
    cd ${run_path}/jobs-dSQ_makeGAM    

    # step .2: submit the dSQ jobs
    sbatch ${run_path}/jobs-dSQ_makeGAM/dsq-${data_prefix}*
    echo "dSQ jobs submitted for ${data_prefix} (run fitGAM in parallel)"
    echo "after the jobs are finished running, submit an additonal job to complete the analysis"
fi 
