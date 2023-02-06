rm(list = ls()); gc()
set.seed(42)

#################################
####### Pre-Configuration #######
#################################

args = commandArgs(trailingOnly = TRUE)
data_prefix = args[1]
run_path = args[2]
data_path = args[3]
use_coding = args[4]

setwd(run_path)
# age cutoff & samples cutoff
lb = 20
sample_keep = "Control"

cat("\n")
cat("---- Run makeSeurat&trajAnalysis for : ", data_prefix, "----\n")

#################################
#################################



cat("\n")
cat("\n")
cat("\n")
cat("##################################################################\n")
cat("######################## Begin makeSeurat ########################\n")
cat("##################################################################\n")
cat("\n")
cat("\n")
cat("\n")


#### Initialization ####

### helper function

source("0_helper_functions.R")


### load packages

library(SingleCellExperiment)
library(SeuratObject)
library(SeuratDisk)
library(Seurat)

library(zellkonverter)
library(reticulate)
library(anndata)
library(SparseM)

library(harmony)
library(stringr)

library(readr)
library(slingshot)
library(viridis)

library(MASS)


### load DATA 

# load seurat object (latest loading method - Jan 12 2023)
seurat_obj_query <- read_h5ad( paste0(data_path, data_prefix, "_annotated.h5ad") )
count.data <- t(seurat_obj_query$raw$X)
count.data <- as(count.data,"matrix.csr")
count.data <- as(count.data,"dgCMatrix")
colnames(count.data) <- row.names(seurat_obj_query$obs)
row.names(count.data) <- row.names(seurat_obj_query$var)
seurat_obj_query <- SeuratObject::CreateSeuratObject(counts = count.data, meta.data = seurat_obj_query$obs)

seurat_obj_query$predictions = seurat_obj_query$azimuth
Idents(seurat_obj_query) <- seurat_obj_query$predictions

# read in only the coding features
coding <- read.delim("./Protein_coding_genes_biomart.txt")


cat("\n")
cat("---- Data loaded for : ", data_prefix, "----\n")
cat("\n")
cat("print initial cell type dimension\n")
print(dim(seurat_obj_query))
cat("\n")
cat("print initial cell type distribution\n")
print(table(seurat_obj_query$predictions))
cat("\n")



#### Subsetting ####

cat("\n")
cat("---- Start subsetting : ", data_prefix, "----\n")
cat("\n")

# 1. subset based on samples 
#    a. only IT neurons
#    b. age >= 20
#    c. only Control samples

# extract only IT neurons 
seurat_obj_query <- subset(seurat_obj_query, 
                           idents = c("L2/3 IT", "L6 IT", "L4 IT", "L5 IT"))

# b. filter by age >= 20 && only Control samples
metadata <- organize_metadata(data_prefix)                          # a helper function
metadata <- match_metadata(data_prefix, metadata, seurat_obj_query) # a helper function

diag_logic = "Diagnosis" %in% names(metadata)
if (diag_logic) {
    keep_age <- which(metadata$Age >= lb & metadata$Diagnosis == sample_keep)
} else {
    keep_age <- which(metadata$Age >= lb)
}

seurat_obj_query$Age = metadata$Age
if (diag_logic) {
  seurat_obj_query$Diagnosis = metadata$Diagnosis
}

seurat_obj_query <- seurat_obj_query[, keep_age]


cat("-- After subset by Age&Control&CellTypes, data dimension: ", dim(seurat_obj_query), "\n")
cat("\n")
cat("Examine & save QC results: ", dim(seurat_obj_query), "\n")
cat("\n")

pdf(file=paste0("./results/", data_prefix, "_QC_after_age_ctrl_celltype.pdf"))
truehist(seurat_obj_query$nCount_RNA)
truehist(seurat_obj_query$nFeature_RNA)
VlnPlot(seurat_obj_query, c("nCount_RNA", "nFeature_RNA"), pt.size = 0)
dev.off()

cat("QC results saved\n")

 
# 2. filter out low quality cells (not set upper limit first)

seurat_obj_query <- subset(x = seurat_obj_query, subset = nCount_RNA > 1000  & nFeature_RNA > 500)


# 3. preserve only coding features (optional)
if (use_coding == TRUE) {
  keep_gene <- which( rownames(seurat_obj_query) %in% coding$hgnc_symbol )
  seurat_obj_query <- seurat_obj_query[keep_gene, ]
}


# this is actually important to do
seurat_obj_query$predictions = as.character(seurat_obj_query$predictions)
seurat_obj_query$individualID = as.character(seurat_obj_query$individualID)

cat("-- After full subsetting, data dimension: ", dim(seurat_obj_query), "\n")
cat("\n")


cat("Cell types distribution after subset:\n")
print(table(seurat_obj_query$predictions))
cat("\n")

cat("Batch distribution after subset: \n")
print(table(seurat_obj_query$individualID))
cat("\n")

cat("Age distribution after subset:\n")
print(table(seurat_obj_query$Age))
cat("\n")

if (diag_logic) {
    cat("Diagnosis distribution after subset:\n")
    print(table(seurat_obj_query$Diagnosis))
    cat("\n")
} else {
    cat("No diagnosis info for", data_prefix, "\n")
    cat("  taking it as only having Control samples\n")
    cat("\n")
}



cat("\n")
cat("---- Finish subsetting : ", data_prefix, "----\n")
cat("\n")






#### Normalization ####

cat("\n")
cat("---- Begin Normalization with SCTransform : ", data_prefix, "----\n")
cat("\n")


# Normalization with SCTransform
seurat_obj_query <- SCTransform(seurat_obj_query, 
                                variable.features.n = NULL, 
                                variable.features.rv.th = 1.3, 
                                conserve.memory=T)

cat("\n")
cat("---- Finish Normalization : ", data_prefix, "----\n")
cat("\n")
cat("Print data dimension after Normalization\n")
print(dim(seurat_obj_query))
cat("\n")





#### Batch effect Correction ####

cat("\n")
cat("---- Begin Batch Effect Correction: ", data_prefix, "----\n")
cat("\n")


seurat_obj_query$samples = seurat_obj_query$individualID

# IMPORTANT step: run PCA before Harmony
seurat_obj_query <- RunPCA(seurat_obj_query)
seurat_obj_query <- RunHarmony(seurat_obj_query, "samples", 
                               plot.convergences=T, project.dim=F,
                               max.iter.harmony = 50)

cat("\n")
cat("---- Finish Batch Effect Correction: ", data_prefix, "----\n")
cat("\n")





#### Generate Embeddings ####

cat("\n")
cat("---- Begin Generating Embeddings: ", data_prefix, "----\n")
cat("\n")

# (optional) Find Clusters Find Neighbors
# another way of clustering - cluster using harmony (not using SCT_snn)
seurat_obj_query <- FindNeighbors(seurat_obj_query, reduction = "harmony", dims = 1:30)
seurat_obj_query <- FindClusters(seurat_obj_query, resolution = 0.5) 


# UMAP embeddings
seurat_obj_query <- RunUMAP(seurat_obj_query, reduction="harmony", dims = 1:30)

pdf(file=paste0("./results/", data_prefix,"_UMAP_savings.pdf"))

DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5, 
        group.by = "seurat_clusters") + NoLegend()
DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5, 
        group.by = "samples") + NoLegend()
DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5, 
        group.by = "predictions") + NoLegend()

dev.off()


cat("\n")
cat("---- Finish saving results of UMAP Embeddings: ", data_prefix, "----\n")
cat("\n")



#### Another Subset: preserve only high variable genes ####


cat("\n")
cat("---- Subset to preserve only variable genes: ", data_prefix, "----\n")
cat("this is to save computation time for the downstream analysis\n")
cat("\n")

seurat_obj_query <- subset(x = seurat_obj_query,
                           features=VariableFeatures(object=seurat_obj_query))
# to save computation time for the following analysis

cat("\n")
cat("---- Finish Subset by variable genes: ", data_prefix, "----\n")
cat("\n")
cat("\n")
cat("---- **** Final data dimension **** ----\n")
cat(dim(seurat_obj_query)[1], "genes\n")
cat(dim(seurat_obj_query)[2], "cells\n")
cat("this is important for parallelly fitting GAM in the following analysis")
cat("\n")




#### Save Processed file ####

cat("\n")
cat("---- Saving Seurat Object ----\n")

save_dir <- paste0("./data/", data_prefix, "_seurat_subset", ".h5seurat")

# in h5Seurat
SaveH5Seurat(seurat_obj_query, 
             filename=save_dir, 
             overwrite = T)

cat("---- Saved in .h5Seurat format ----\n")
cat("\n")

# convert to AnnData
Convert(save_dir, dest = "h5ad")

cat("---- Saved in .h5ad format ----\n")
cat("\n")








cat("\n")
cat("\n")
cat("##################################################################\n")
cat("####################### Begin trajAnalysis #######################\n")
cat("##################################################################\n")
cat("\n")
cat("\n")




cat("\n")
cat("---- Begin trajectory analysis for", data_prefix, "----\n")
cat("dimension of data: ")
cat(dim(seurat_obj_query)[1], "genes\n")
cat(dim(seurat_obj_query)[2], "cells\n")


#### Sling Traj (embedded on UMAP) ####

# extract UMAP
umap <- Embeddings(seurat_obj_query, reduction='umap')

# check seurat Idents - Turn to cell types
Idents(seurat_obj_query) <- seurat_obj_query$predictions


cat("---- Slingshot BEGIN ----\n")
# Run slingshot - Construct Trajectories
seurat_sling <- slingshot(as.matrix(umap), 
                          clusterLabels = Idents(seurat_obj_query), 
                          start.clus='L2/3 IT')
cat("---- Slingshot DONE  ----\n")



# visualize trajectories
plotcol <- viridis::turbo(length(levels(seurat_obj_query)))
col_vec <- rep("", ncol(seurat_obj_query))
col_vec[as.character(Idents(seurat_obj_query)) == "L2/3 IT"] <- plotcol[1]
col_vec[as.character(Idents(seurat_obj_query)) == "L4 IT"]   <- plotcol[2]
col_vec[as.character(Idents(seurat_obj_query)) == "L5 IT"]   <- plotcol[3]
col_vec[as.character(Idents(seurat_obj_query)) == "L6 IT"]   <- plotcol[4]

pdf(file=paste0("./results/", data_prefix, "_UMAP_traj.pdf"))
plot(umap[,1], umap[,2], col=col_vec, pch = 16, cex = 0.6)
legend("bottomright",
       legend=c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT" ), 
       col=plotcol, pch=16, bty='n', cex=0.6)
lines(SlingshotDataSet(seurat_sling), lwd=2, col='black')
dev.off()


saveRDS(seurat_sling, paste0("./data/", data_prefix, 
                             "_seurat_subset_sling", ".rds"))

cat("\n")
cat("---- trajAnalysis results SAVED ----\n")




