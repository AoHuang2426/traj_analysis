## Pipeline Manual: Trajectory Analysis across Cortical Layers

### Run the Pipeline in One Command

To run the pipeline, configure the script `run_analysis.sh` first. Things to configure:

+   `data_prefix`: the data prefix of the dataset (e.g. Kriegstein)
+   `run_path`: the present working directory that contains these codes
+   `data_path`: the directory to the datasets that are required in this analysis. For now it is the `Hybrid_cell_scheme` directory.
+   `use_coding`: boolean value to decide whether or not only the coding features are used in the analysis. Default is `FALSE` meaning that all features are included. 
+   `run_parallel`: boolean value to decide whether or not to run fitGAM in parallel. If we want to complete the analysis in a single command, we need to set it to `FALSE`, but it might take longer time. 

After configuring these parameters, simply submit the job to Slurm, which will complete the trajectory analysis and fitGAM analysis. The Seurat objects are GAM objects will be saved in `./data`, and all the results will be saved in `./results`.

```shell
sbatch ./run_analysis.sh
```

After the pipeline is finished, run the `4_postprocessing_dsanalysis.R` interactively to produce the smoothly varying genes, heatmap visualization and common genes across datasets. Note that we will need to specify which datasets are to be included.

### Run the fitGAM in Parallel (optional)

Some datasets are large with many cells and features (e.g. Ma_Sestan, UCLA-ASD). Set `run_parallel=TRUE` would utilize the `dSQ` job arrays system to run fitGAM in parallel. Doing this requires additional steps but might save some running time. 

#### Step 1

Set `run_parallel=TRUE` in `run_analysis.sh` and submit the job to Slurm.

```
sbatch ./run_analysis.sh
```

#### Step 2

When the `dSQ` jobs are finished running, submit an addtional script `after_fitGAM_parallel.sh` to Slurm to complete the pipeline. Make sure to keep the `data_prefix` and `run_path` consistent as in the `run_analysis.sh` script.

```
sbatch ./after_fitGAM_parallel.sh
```

#### Step 3

After the pipeline is finished, run `4_postprocessing_dsanalysis.R` interactively to produce and visualize the results. 



**Credits** to Tiernon Riesenmy ([link](https://github.com/TiernonRR/Trajectory_Analysis)), Dr. Cagatay Dursun and Dr. Prashant Emani.
