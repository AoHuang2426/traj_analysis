# output metadata by the order of ID | Sex | Age | PMI | Diagnosis (optional)
organize_metadata <- function (data_prefix) {
  
  data_path = "/gpfs/gibbs/pi/gerstein/pse5/Hybrid_cell_scheme/"
  
  # read in reduced meta data
  if (data_prefix == "IsoHuB") {
    metadata <- read.csv(paste0(data_path, data_prefix, "_metadata_reduced.csv"))
  } else if (data_prefix == "CMC") {
    metadata <- read.delim(paste0(data_path, data_prefix, "_reduced_individual_metadata.tsv"))
  } else {     # "Kriegstein", "Ma_Sestan", "DevBrain", "UCLA-ASD", "SZBD-Kellis"
    metadata <- read.delim(paste0(data_path, data_prefix, "_reduced_metadata.tsv"))
  }
  
  # organize
  names(metadata)[names(metadata) == "individualID" | names(metadata) == "Sample.Library.No."] = "ID"
  names(metadata)[names(metadata) == "Gender" | names(metadata) == "Sex" | names(metadata) == "reportedGender"] = "Sex" 
  names(metadata)[names(metadata) == "Age" | names(metadata) == "ageOfDeath" | names(metadata) == "ageDeath"] = "Age" 
  diag_logic = "Diagnosis" %in% names(metadata) | "primaryDiagnosis" %in% names(metadata)
  
  if (diag_logic) {
    names(metadata)[names(metadata) == "Diagnosis" | names(metadata) == "primaryDiagnosis"] = "Diagnosis"
    meta_output = data.frame(ID = as.character(metadata$ID),
                             Sex = as.character(metadata$Sex),
                             Age = metadata$Age,
                             PMI = metadata$PMI,
                             Diagnosis = as.character(metadata$Diagnosis))
  } else {
    meta_output = data.frame(ID = as.character(metadata$ID),
                             Sex = as.character(metadata$Sex),
                             Age = metadata$Age,
                             PMI = metadata$PMI)
  }
  
  if (data_prefix == "CMC") {
    meta_output$Age = as.numeric(substring(metadata$Age, 1, 2))
  }
  
  # specify the "Control" entries
  if (diag_logic) {
    ind = which(meta_output$Diagnosis == "Not Applicable" | 
                  meta_output$Diagnosis == "control" | 
                  meta_output$Diagnosis == "Control")
    
    meta_output$Diagnosis[ind] = "Control"
  }
  
  return(meta_output)
}


# assign each cell ID with Sex, Age, PMI
match_metadata <- function (data_prefix, metadata, seurat) {
  
  donors <- as.character(seurat$individualID) # for all 4+3 datasets this is equivalent to samples
  
  # special treatment to Ma_Sestan and DevBrain
  if (data_prefix == "Ma_Sestan" | data_prefix == "DevBrain") {
    donors <- sapply(donors, function(x){unlist(strsplit(x, "_"))[1]})
  }
  
  diag_logic = "Diagnosis" %in% names(metadata)
  
  if (diag_logic) {
    meta_output <- data.frame(matrix(0, nrow = ncol(seurat), ncol = 4)) # col: Sex, Age, PMI, Diagnosis
    # row: # of cells in seurat
    names(meta_output) = c("Sex", "Age", "PMI", "Diagnosis")
    meta_output$Diagnosis = as.character(meta_output$Diagnosis)
    rownames(meta_output) <- colnames(seurat)
    
    # nrow(metadata) is unique(donors)
    for (i in 1:nrow(metadata)) {
      ind <- which(donors == metadata$ID[i])
      if (length(ind) > 0) { # include only the donor id existed in the seurat obj
        meta_output[ind, ] = metadata[i, c("Sex", "Age", "PMI", "Diagnosis")]
      }
    }
    
  } else {
    meta_output <- data.frame(matrix(0, nrow = ncol(seurat), ncol = 3)) # col: Sex, Age, PMI
    # row: # of cells in seurat
    names(meta_output) = c("Sex", "Age", "PMI")
    rownames(meta_output) <- colnames(seurat)
    
    # nrow(metadata) is unique(donors)
    for (i in 1:nrow(metadata)) {
      ind <- which(donors == metadata$ID[i])
      if (length(ind) > 0) { # include only the donor id existed in the seurat obj
        meta_output[ind, ] = metadata[i, c("Sex", "Age", "PMI")]
      }
    }
  }
  
  return(meta_output)
  
}