rm(list = ls());gc()

#################################
####### Pre-Configuration #######
#################################

setwd("~/project/proj_spatialCon_rt2/traj_analysis/")
datasets = c("Kriegstein")

library(reticulate)
library(slingshot)
library(tradeSeq)
library(SeuratDisk)
library(Seurat)
library(reticulate)
library(SingleCellExperiment)
library(pheatmap)
library(UpSetR)
library(stringr)
library(cowplot)

#################################
#################################




for (data_prefix in datasets) {
  
  ### helper function
  heatMap <- function(lineage, title="continuousGeneExpr", rows=F){
    yhatSmooth <- predictSmooth(gam, gene=lineage, nPoints=50, tidy=F)
    heatSmooth <- pheatmap( t(scale(t(yhatSmooth[, 1:50]))), 
                            cluster_cols=F, show_rownames=rows, show_colnames=F,
                            # filename = paste0("./results/",data_prefix, "_" , title, ".png") 
                            )
  } 
  
  
  ### load data
  seurat_sling <- readRDS(paste0("./data/", data_prefix, 
                                 "_seurat_subset_sling", ".rds"))
  gam <- readRDS(paste0("./data/gam/", data_prefix, "_gam/",
                        data_prefix, "_gam_it_all.rds"))
  
  
  #### remove high expr genes in ONLY one cluster ####
  
  to_remove <- list()
  levels <- c("L2_3_IT", "L4_IT", "L5_IT", "L6_IT")
  for (clust in levels){
    if(data_prefix == "Kriegstein"){
      t <- read.table(file=paste0("/home/tr455/project/DEGs_trajectory_Analysis/Kriegstein_", 
                                  clust, "_DEGs.tsv"), sep="\t", header=T)
    } else if (data_prefix == "CMC" | data_prefix == "SZBD-Kellis" | data_prefix == "UCLA-ASD" ) { 
      t <- read.table(file=paste0("/gpfs/gibbs/pi/gerstein/pse5/Hybrid_cell_scheme/",
                                  data_prefix, "_", clust, "_DEGs.csv"), 
                      sep="\t", header=T)
    } else {
      t <- read.table(file=paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
                                  data_prefix, "_", clust, "_DEGs.csv"), 
                      sep="\t", header=T)
    }
    indices <- row.names(t[which(t['p_val_adj'] < 0.05),])
    to_remove <- append(to_remove, indices)
  }
  
  to_remove <- unique(to_remove)
  sub <- rownames(gam) %in% to_remove
  gam <- gam[which(sub==F),]
  cat("----", dim(gam)[1], "genes remaining for", data_prefix, "----\n")

  
  
  #### Association Test ####
  
  rowData(gam)$assocRes <- associationTest(gam, l2fc = log2(2), lineages=T)
  assocRes <- rowData(gam)$assocRes
  
  if (data_prefix == "Kriegstein" | data_prefix == "CMC" | data_prefix == "SZBD-Kellis" | data_prefix == "IsoHuB") { 
    
    l = rownames(assocRes)[ which( p.adjust(assocRes$pvalue_1, "fdr") <= .05) ]
    
  } else if (data_prefix == "Ma_Sestan") {
    # conds = FALSE
    l1 <- rownames(assocRes)[ which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05) ] 
    
  } else if (data_prefix == "DevBrain") {
    # conds = TRUE
    l1_cond <- rownames(assocRes)[ which(p.adjust(assocRes$`pvalue_lineage1_conditionNot Applicable`, "fdr") <= 0.05) ]
  }
  
  
  
  
  #### Remove NAs, Infs and Zeros expr before Plot ####
  
  if (data_prefix == "Kriegstein" | data_prefix == "CMC" | data_prefix == "SZBD-Kellis" | data_prefix == "IsoHuB") { 
    to_check <- l
  } else if (data_prefix == "Ma_Sestan") {
    to_check <- l
  }
  
  
  # return a gene x value matrix, where dim(value) = 100 = 50(lineage1) + 50(lineage2)
  yhatSmooth <- predictSmooth(gam, gene=to_check, nPoints=50, tidy=F)
  scale_result <- t( scale(t(yhatSmooth[, 1:50])) ) # scale for lineage 1
  
  infs <- rep("", length(rownames(scale_result))) # num of genes 
  NAs <- rep("", length(rownames(scale_result)))
  zeros <- rep("", length(rownames(scale_result)))
  
  # check 0 expr genes
  for(i in 1:length(to_check)){
    gene <- to_check[i]
    if( all(yhatSmooth[gene, ] == 0) ){
      zeros[i] <- gene
    }
  }
  
  # check NA & infs expr genes
  for(i in 1:length(to_check)){
    gene <- to_check[i]
    
    if(any(is.na(scale_result[gene,]) )){
      NAs <- gene
    } else if(any(scale_result[gene,] == Inf)){
      infs[i] <- gene
    }
  }
  
  cat("-- There are", unique(infs),  "inf expr --\n")
  cat("-- There are", unique(NAs),   "NAs expr --\n")
  cat("-- There are", unique(zeros), "zeros expr --\n")
  
  # remove NA, Inf, zeros
  to_remove <- union(unique(infs), union(unique(NAs), unique(zeros)))
  to_remove <- unique(to_remove)
  sub <- rownames(gam) %in% to_remove
  gam <- gam[which(sub==F),]
  cat("----", dim(gam)[1], "genes remaining for", data_prefix, "(after filter out Infs NAs & Zeros) ----\n")
  # Kriesgstein: 1345 genes remaining
  # IsoHuB:      1526 genes remaining
  # CMC:         895 genes remaining
  # SZBD-Kellis: 994 genes remaining 
  
  
  
  #### Re-do Association Test (after removal) ####
  
  rowData(gam)$assocRes <- associationTest(gam, l2fc = log2(2), lineages=T)
  assocRes <- rowData(gam)$assocRes
  
  if (data_prefix == "Kriegstein" | data_prefix == "CMC" | data_prefix == "SZBD-Kellis" | data_prefix == "IsoHuB") { 
    
    l = rownames(assocRes)[ which( p.adjust(assocRes$pvalue_1, "fdr") <= .05) ]
    
  } else if (data_prefix == "Ma_Sestan") {
    # conds = FALSE
    l1 <- rownames(assocRes)[ which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05) ] 
    
  } else if (data_prefix == "DevBrain") {
    # conds = TRUE
    l1_cond <- rownames(assocRes)[ which(p.adjust(assocRes$`pvalue_lineage1_conditionNot Applicable`, "fdr") <= 0.05) ]
  }
  
  
  
  
  #### output ####
  
  if (data_prefix == "Kriegstein" | data_prefix == "CMC" | data_prefix == "SZBD-Kellis" | data_prefix == "IsoHuB") { 
    heatMap(lineage = l)
    saveRDS(l, paste0("./data/ds_analysis/", data_prefix, "_l.rds"))
    
  } else if (data_prefix == "Ma_Sestan") {
    heatMap(l1, "l1")
    # saveRDS(l1, paste0("./data/ds_analysis/", data_prefix, "_l1"))
  }
  
  
}














#### Common genes analysis ####

l1_k <- readRDS("./data/ds_analysis/Kriegstein_l.rds")
l1_i <- readRDS("./data/ds_analysis/IsoHuB_l.rds")
l1_c <- readRDS("./data/ds_analysis/CMC_l.rds")
l1_s <- readRDS("./data/ds_analysis/SZBD-Kellis_l.rds")

u_k <- unique(l1_k)
u_i <- unique(l1_i)
u_c <- unique(l1_c)
u_s <- unique(l1_s)

common_genes <- intersect(u_k, intersect(u_i, intersect(u_c, u_s)))

pdf(file=paste0("./results/commonGenesVenn.pdf"))
UpSetR::upset(fromList(list(Kriegstein=u_k, IsoHuB=u_i, CMC=u_c, SZBD = u_s)))
dev.off()


#### GAM with common genes for all datasets ####

gam_K <- readRDS(paste0("./data/gam/", "Kriegstein", "_gam/", "Kriegstein", "1_gam_it_all.rds"))
gam_I <- readRDS(paste0("./data/gam/", "IsoHuB", "_gam/", "IsoHuB", "1_gam_it_all.rds"))
gam_C <- readRDS(paste0("./data/gam/", "CMC", "_gam/", "CMC", "1_gam_it_all.rds"))
gam_S <- readRDS(paste0("./data/gam/", "SZBD-Kellis", "_gam/", "SZBD-Kellis", "1_gam_it_all.rds"))

### Kriegstein
yhatSmooth <- predictSmooth(gam_K, gene=common_genes, nPoints=50, tidy=F)
heatSmooth_k <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, cluster_rows = T, clustering_distance_rows = 'euclidean', clustering_method = 'ward.D2',
                       show_rownames=T,show_colnames=F, fontsize=6,
                       filename = paste0("./results/", "Kriegstein", "_" , "common_genes", ".png")) 

# order of genes in Kriegstein
row_gene_k = heatSmooth_k[["tree_row"]][["labels"]]
row_gene_k = row_gene_k[heatSmooth_k[["tree_row"]][["order"]]]


### IsoHuB
yhatSmooth <- predictSmooth(gam_I, gene=common_genes, nPoints=50, tidy=F)
heatSmooth_i <- pheatmap(t(scale(t(yhatSmooth[row_gene_k, 1:50]))), cluster_cols=F, cluster_rows = F,
                       show_rownames=T,show_colnames=F, fontsize=6,
                       filename = paste0("./results/", "IsoHuB", "_" , "common_genes", ".png"))

### SZBD
yhatSmooth <- predictSmooth(gam_S, gene=common_genes, nPoints=50, tidy=F)
heatSmooth_i <- pheatmap(t(scale(t(yhatSmooth[row_gene_k, 1:50]))), cluster_cols=F, cluster_rows = F,
                         show_rownames=T,show_colnames=F, fontsize=6,
                         filename = paste0("./results/", "SZBD-Kellis", "_" , "common_genes", ".png"))

### CMC
yhatSmooth <- predictSmooth(gam_C, gene=common_genes, nPoints=50, tidy=F)
heatSmooth_c <- pheatmap(t(scale(t(yhatSmooth[row_gene_k, 1:50]))), cluster_cols=F, cluster_rows = F,
                       show_rownames=T,show_colnames=F, fontsize=6,
                       filename = paste0("./results/", "CMC", "_" , "common_genes", ".png"))



# # heatmap together
# yhatSmooth_k <- predictSmooth(gam_K, gene=common_genes, nPoints=50, tidy=F)
# yhatSmooth_i <- predictSmooth(gam_I, gene=common_genes, nPoints=50, tidy=F)
# yhatSmooth_c <- predictSmooth(gam_C, gene=common_genes, nPoints=50, tidy=F)
# 
# matrix_all = merge(t(scale(t(yhatSmooth_k[, 1:50]))),
#                    t(scale(t(yhatSmooth_i[, 1:50]))),
#                    t(scale(t(yhatSmooth_c[, 1:50]))), by = "row.names", all = T)


### save common genes
saveRDS(common_genes, paste0("./data/ds_analysis/", "common_genes.rds"))


# fileConn <- file("/home/tr455/project/coding/common_genes_it.txt")
# writeLines(common_genes, fileConn)
# close(fileConn)
#
# gene_union <- c(unique(rownames(gam_D)),unique(rownames(gam_M)),
#   unique(rownames(gam_K)), unique(rownames(gam_I)))
#
# fileConn <- file("/home/tr455/project/coding/union_genes_it.txt")
# writeLines(gene_union, fileConn)
# close(fileConn)

### (optional) same row names, 4 plots together




#### GO Enrichment analysis ####

# load package

# BiocManager::install("clusterProfiler")  #用来做富集分析
# BiocManager::install("topGO")  #画GO图用的
# BiocManager::install("Rgraphviz")
# # BiocManager::install("pathview") #看KEGG pathway的
# BiocManager::install("org.Hs.eg.db") #这个包里存有人的注释文件

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
# library(pathview)
library(org.Hs.eg.db)

# load data

common_genes <- readRDS(paste0("./data/ds_analysis/", "common_genes.rds"))
# common_genes_t <- read.table("common_genes_it_bytiernon.txt", header=FALSE)
# common_genes_t <- as.character(common_genes_t$V1)

# # GO analysis
# 
# DEG.gene_symbol = common_genes
# DEG.entrez_id = mapIds(x = org.Hs.eg.db,
#                        keys = DEG.gene_symbol,
#                        keytype = "SYMBOL",
#                        column = "ENTREZID")
# DEG.entrez_id = na.omit(DEG.entrez_id)
#   # tiernon: 360 -> 263
#   # Ao: 100 -> 100
#   # tiernon: 10% non-coding features
# 
# erich.go.BP = enrichGO(gene = DEG.entrez_id,
#                        OrgDb = org.Hs.eg.db,
#                        keyType = "ENTREZID",
#                        ont = "MF",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.05)
# 
# dotplot(erich.go.BP)
# 
# common_genes <- data.frame(common_genes)
# 
# 
# 
# 
# 
# 
# 
# 
