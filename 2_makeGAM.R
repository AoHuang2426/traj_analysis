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
cat("######################### Begin makeGAM ##########################\n")
cat("##################################################################\n")
cat("\n")
cat("\n")
cat("\n")


#### helper function

source("0_helper_functions.R")


#### Load packages #### 

library(SeuratDisk)
library(Seurat)
library(reticulate)
library(tradeSeq)
library(SingleCellExperiment, quietly = T)
library(stringr)
library(BiocParallel)
library(ggplot2)

# read data & DO SUBSETTING

cat("\n")
cat("---- Reading in data for", data_prefix, "----\n")
cat("\n")

filepath <- paste0("./data/", data_prefix, "_seurat_subset", ".h5seurat")
seurat_obj_query <- LoadH5Seurat(filepath)
Idents(seurat_obj_query) <- seurat_obj_query$predictions

cat("-- Data dimension: ", dim(seurat_obj_query), "\n")
cat("\n")

cat("-- Cell types distribution: \n")
print(table(Idents(seurat_obj_query)))


#### metadata info check ####

# NOTE: not every sample info in metadata is in seurat$samples --> changed to seurat$individualID
metadata <- organize_metadata(data_prefix)

# check donor (sample) info
sample_counts = table(seurat_obj_query$individualID)

pdf(file=paste0("./results/", data_prefix, "_sample_histogram.pdf"))
plot(sort(sample_counts, decreasing = T), 
     type="h", las=2, cex.axis=0.60, 
     ylab="Cell Counts")
dev.off()

cat("---- Data Samples/Individuals distribution SAVED ----\n")
cat("\n")


cat("\n")
cat("---- Constructing Design Matrix ----\n")
cat("\n")
cat("-- Covariates: Sex, Age, PMI, Props, dcounts \n")
cat("Props:   the cell type fraction in the sample that a cell belongs to \n")
cat("dcounts: the sample size (cell number) of the sample that a cell belongs to \n")
cat("\n")

#### Design matrix of mixed effects ####

# Sex, Age, PMI, sample counts, sample fractions

# sample's contribution (fraction) of each cell types
sample_type <- unique(seurat_obj_query$individualID)
frac_per_cell_type <- sapply(sample_type, function(x){
  table( Idents(seurat_obj_query)[which(seurat_obj_query$individualID == x)] )})
frac_per_cell_type <- t( frac_per_cell_type )
frac_per_cell_type <- prop.table(frac_per_cell_type, 2)
frac_per_cell_type <- data.frame(frac_per_cell_type)
frac_per_cell_type$samples <- rownames(frac_per_cell_type)

# new index of cell types (to match that with frac_per_cell_type)
ids <- sapply(as.character(Idents(seurat_obj_query)), function(x) {str_replace_all(x, " ", ".")} )
ids <- sapply(ids, function(x) str_replace_all(x, "/", "."))

# match samples with cell types 
tmp <- cbind(seurat_obj_query$individualID, ids)

# construct design matrix
U <- match_metadata(data_prefix, metadata, seurat_obj_query)
U$dcounts <- sapply(seurat_obj_query$individualID, function(x){sample_counts[x]}) # the cell number of this sample that this cell belong to (the cell type)
lst <- rep(0, length(ids))
for(i in 1:length(ids)){
  lst[i] <- frac_per_cell_type[ tmp[i,][[1]], tmp[i,][[2]] ]
}
U$props <- lst                                                                    # the cell fraction of this cell type in this sample that this cell belong to (the cell type in this sample)

# binarize Sex input
U$Sex <- lapply(U$Sex, function(x) {as.integer(tolower(x) == "male")})
U$Age <- lapply(U$Age, function(x) {as.integer(x)})
U$PMI <- lapply(U$PMI, function(x) {as.double(x)})


cat("\n")
cat("---- Configure Model ----\n")
cat("\n")

#### Model ####

sex <- as.factor(unlist(U$Sex)) 
pmi <- unlist(U$PMI)
pmi[is.na(pmi)] <- mean(pmi[!is.na(pmi)])
age <- unlist(U$Age)
prop <- unlist(U$props)
dcount <- unlist(U$dcounts)

options(na.action="na.pass")
if (data_prefix == "DevBrain") {
  U_final <- model.matrix(~ age + pmi + prop + dcount)  # only has male, no female       
} else {
  U_final <- model.matrix(~ age + pmi + sex + prop + dcount)  
}
     

counts <- seurat_obj_query@assays$SCT@counts



#### save for next step #### 
saveRDS(U_final, paste0("./data/gam/", data_prefix, "_U_it.rds"))
saveRDS(counts,file=paste0("./data/gam/", data_prefix, "_counts_it.rds"))


cat("---- Design Matrix saved for", data_prefix, "----\n")
cat("\n")  
  