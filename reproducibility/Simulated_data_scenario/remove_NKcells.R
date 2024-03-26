# This file is used to remove NK cells from the artificially-mixed tissue-blood data

library(dplyr)
library(Seurat)
library(patchwork)

# Path to the input file
patht_to_artificially_mixed_TB <- "E:/Tays_workspace/data/data_merged_tissue_blood_residents_DEMO_noNKcellTissue.rds"

data.artificially.mixed.TB <- readRDS(patht_to_artificially_mixed_TB)

# Filter out NK cells
list_approved_cell_types <- c("Epithelial cell", "T-cell", "Monocyte", "B-cell", "Tissue stem cell")

data.artificially.mixed.TB.filtered <- subset(x = data.artificially.mixed.TB, subset = final_annotation %in% list_approved_cell_types)

# Save data
saveRDS(data.artificially.mixed.TB.filtered, "E:/Tays_workspace/data/data_merged_tissue_blood_residents_DEMO_noNKcellTissueBlood.rds")