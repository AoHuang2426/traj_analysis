rm(list = ls())
gc()


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
cat("############### Begin combineGAM (parallel fitGAM) ###############\n")
cat("##################################################################\n")
cat("\n")
cat("\n")
cat("\n")


library(SingleCellExperiment)
library(tradeSeq)

# get gene numbers
counts <- readRDS(paste0("./data/gam/", data_prefix, "_counts_it.rds"))
n_gene <- dim(counts)[1]

cat("---- Start combining", data_prefix, "parallel gam files ----\n")
cat("\n")

# load in the FIRST file

path <- paste0("./data/gam/", data_prefix, "_gam/")
combined_GAM <- readRDS(paste0(path, data_prefix, "1_gam_it.rds"))
# load in the REMAINING file

chunk_length <- 200
to_split <- (chunk_length+1):ngenes
chunks <- split(to_split, ceiling(seq_along(to_split)/chunk_length))

for(chunk in chunks) {
  c_min <- min(chunk)
  gam <- readRDS( paste0(path, data_prefix, c_min, "_gam_it.rds") )
  combined_GAM <- rbind(combined_GAM, gam)
}


saveRDS(combined_GAM, paste0("./data/gam/", data_prefix, "_gam/",
                           data_prefix, "_gam_it_all.rds"))

cat("---- saved", data_prefix, "----\n\n")
