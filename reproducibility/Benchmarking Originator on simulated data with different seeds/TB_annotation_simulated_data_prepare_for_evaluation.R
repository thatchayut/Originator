# This file is used to prepare Originator results on simulated data for evaluation
#   using different metrics in python

library(Seurat)
library(dplyr)

setwd("~/temp/soft_link_nfs_originator_GB_revision")

# Visualization
# seeds <- c(0)
# celltypes <- c("tcell")
# unified_celltypes <- c("T-cell")

seeds <- c(0, 42, 64, 123, 894)
celltypes <- c("tcell", "mono", "b")
unified_celltypes <- c("T-cell", "Monocyte", "B-cell")

for (i in 1:length(celltypes)) {
  
  celltype <- celltypes[[i]]
  unified_celltype <- unified_celltypes[[i]]
  
  message(celltype)
  
  for (seed in seeds) {
    
    data <- readRDS(file.path("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/simulatedData_seeds", celltype, seed, paste0(unified_celltype, "_annotated.rds")))
    
    # Extract to get only Query
    data_query <- subset(data, subset = orig.ident == "Query")
    
    # Extract required information
    data_query_meta <- data_query@meta.data
    
    # Prepare data in the required format for the evaluation code in python
    data_for_evaluation <- data_query_meta %>% select(orig.ident, final_annotation,
                                                      origin_tb, origin_groundtruth)
    
    # Add column for binary class value
    data_for_evaluation$TB_origin_class <- "NA"
    data_for_evaluation$TB_origin_class[data_for_evaluation$origin_tb == "Tissue"] <- 1
    data_for_evaluation$TB_origin_class[data_for_evaluation$origin_tb == "Blood"] <- 0
    
    data_for_evaluation$origin_groundtruth_class <- "NA"
    data_for_evaluation$origin_groundtruth_class[data_for_evaluation$origin_groundtruth == "Tissue"] <- 1
    data_for_evaluation$origin_groundtruth_class[data_for_evaluation$origin_groundtruth == "Blood"] <- 0
    
    # Change column name 
    colnames(data_for_evaluation) <- c("orig.ident", "final_annotation", 
                                       "TB_origin", "origin_groundtruth",
                                       "TB_origin_class", "origin_groundtruth_class")
    
    # Save output file
    write.csv(data_for_evaluation, file.path("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/simulatedData_seeds", celltype, seed, paste0("prediction_TB_", unified_celltype, "_TBannotation.csv")))
    
  }
}



