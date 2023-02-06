rm(list = ls()); gc()

#################################
####### Pre-Configuration #######
#################################

args = commandArgs(trailingOnly = TRUE)
data_prefix = args[1]
run_path = args[2]

setwd(run_path)

cat("-- data_prefix: ", data_prefix, "\n")
cat("-- runn_path:   ", run_path, "\n")

#################################
#################################

cat("\n")
cat("\n")
cat("\n")
cat("##################################################################\n")
cat("################ Begin makeJob (parallel fitGAM) #################\n")
cat("##################################################################\n")
cat("\n")
cat("\n")
cat("\n")



cat("---- make GAM job for ", data_prefix, "----\n")
cat("\n")
  
# get gene numbers
counts <- readRDS(paste0(run_path, "/data/gam/", data_prefix, "_counts_it.rds"))
cat("-- data dimension: ", dim(counts), "\n")
n_gene <- dim(counts)[1]
n_gene <- as.numeric(n_gene)
cat("\n-- number of variable genes: ", n_gene, "\n")
cat(is.numeric(n_gene), "\n\n")
  
to_split <- 1:n_gene
chunk_length <- 200
chunks <- split( to_split, ceiling(seq_along(to_split)/chunk_length) )
  
vect <- rep("", length(chunks))
index <- 1

for(chunk in chunks) {
  c_min <- min(chunk)
  c_max <- max(chunk)
  vect[index] <- paste("module load miniconda; conda init bash; source activate /gpfs/ysm/project/gerstein/ah2426/conda_envs/R4.2.0; Rscript", 
                       paste0(run_path, "/parallel_run/" ,"2.1_fitGAM_parallel.R"), 
                       data_prefix, c_min, c_max, run_path) # the 4 args to be inputed in 2.1_fitGAM_parallel.R       
  index <- index + 1
}
  
filepath <- file(paste0("./data/gam/", data_prefix, "_joblist_parallel.txt"))
writeLines(vect, filepath)

cat("---- makeGAM joblist made for", data_prefix, "----\n")
cat("\n")

close(filepath)



