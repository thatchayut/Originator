# This file is used to generate UMAP plot for the recovery cell types from the artifially mixed blood-tissue resident cell

library(Seurat)

setwd("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref")

'%!in%' <- function(x,y)!('%in%'(x,y))

# Read input data
data <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood.rds")

data$TB_origin_wholeblood <- factor(data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/UMAP_Orignator_TB_recovery_from_simulatedData.pdf", width = 12, height = 10)
DimPlot(data, group.by = "TB_origin_wholeblood")
dev.off()