## Pipeline Manual: Trajectory Analysis across Cortical Layers

### Run the Pipeline in a Single Command

To run the pipeline, configure the script `run_analysis.sh` first. Things to configure:

+   `data_prefix`: the data prefix of the dataset (e.g. Kriegstein)
+   `run_path`: the present working directory that contains these codes
+   `data_path`: the directory to the datasets that are required in this analysis. For now it is the `Hybrid_cell_scheme` directory.
+   `use_coding`: boolean value to decide whether or not only the coding features are used in the analysis. Default is `FALSE` meaning that all features are included. 

After configuring these parameters, simply submit the job to Slurm, which will complete the trajectory analysis and fitGAM analysis. The Seurat objects are GAM objects will be saved in `./data`, and all the results will be saved in `./results`.

```shell
sbatch ./run_analysis.sh
```

After the pipeline is finished, run the `4_postprocessing_dsanalysis.R` interactively to produce the smoothly varying genes, heatmap visualization and common genes across datasets. Note that you will need to specify which datasets are to be inputed.



**Credits** to Tiernon Riesenmy ([link](https://github.com/TiernonRR/Trajectory_Analysis)), Dr. Cagatay Dursun and Dr. Prashant Emani.
