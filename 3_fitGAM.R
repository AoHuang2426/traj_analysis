rm(list = ls()); gc()
set.seed(3)

#################################
####### Pre-Configuration #######
#################################

args = commandArgs(trailingOnly = TRUE)
data_prefix = args[1]
run_path = args[2]

setwd(run_path)

#################################
#################################

cat("\n")
cat("\n")
cat("\n")
cat("##################################################################\n")
cat("########################## Begin fitGAM ##########################\n")
cat("##################################################################\n")
cat("\n")
cat("\n")
cat("\n")


library(tradeSeq)
library(SingleCellExperiment)
library(reticulate)
library(BiocParallel)


cat("\n")
cat("---- Load traj data & design matrix for", data_prefix, "----\n")
cat("\n")

# load data
seurat_sling <- readRDS(paste0("./data/", data_prefix, "_seurat_subset_sling", ".rds"))
U_final <- readRDS(paste0("./data/gam/", data_prefix, "_U_it.rds"))
counts <- readRDS( paste0("./data/gam/", data_prefix, "_counts_it.rds"))
n_gene <- dim(counts)[1]

cat("Design matrix dimension:", dim(U_final) ,"\n")
cat("sct Count matrix dimension:", dim(counts) ,"\n")


cat("\n")
cat("---- fitGAM BEGIN! for", data_prefix, "(not in parallel) ----\n")
cat("\n")

gam <- fitGAM(as.matrix(counts), 
              U=U_final, 
              sds=seurat_sling,
	            genes=1:n_gene) 

cat("\n")
cat("---- fitGAM finished for", data_prefix, " (not in parallel) ----\n")
cat("\n")

saveRDS(gam, file=paste0("./data/gam/", data_prefix, "_gam/" , data_prefix, "_gam_it_all.rds"))



cat("\n")
cat("---- fitGAM saved for", data_prefix, " (all features) ----\n")
cat("\n")



