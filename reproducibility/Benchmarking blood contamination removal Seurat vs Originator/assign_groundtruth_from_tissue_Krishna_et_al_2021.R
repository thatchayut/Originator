# This file is used after "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/derive_groundtruth_from_tissue_Krishna_et_al_2021.R"
# The purpose of this file is to use label transfer results to assign ground truth
#   in ccRCC tissue

# Dataset in this study: Krishna et al. (2021) - Kidney
#   - Use a patient with complete response as an example
#   - Path to file: /nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds
# Given that each patient has paired PBMC and ccRCC tissues
#   1) Use label transfer to identify PBMC in ccRCC tissues by using
#       PBMC as "reference" and ccRCC tissue as "query"
#   2) Treat identified PBMC in ccRCC (+- PBMC) as "blood immune cells" and
#       treat the rest immune cells in CLEANED ccRCC as "tissue-resident immume cells"
#   3) Use these labels as GROUND TRUTH for Originator

##### NOTE: Use Seurat version 4.3.0

library(Seurat)

setwd("~/temp/soft_link_nfs_originator_GB_revision")

# Specify path to output folder
base_output_path <- "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/"

# Read input file
ccRCC <- readRDS("/nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds")

# Read prediction files performed via label transfer in "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/derive_groundtruth_from_tissue_Krishna_et_al_2021.R"
CD4T_predictions <- readRDS(paste0(base_output_path, "label_transfer_result_CD4T.rds"))

CD8T_predictions <- readRDS(paste0(base_output_path, "label_transfer_result_CD8T.rds"))

monocyte_predictions <- readRDS(paste0(base_output_path, "label_transfer_result_Monocyte.rds"))

NK_predictions <- readRDS(paste0(base_output_path, "label_transfer_result_NK.rds"))

Bcell_predictions <- readRDS(paste0(base_output_path, "label_transfer_result_Bcell.rds"))

# Preprocess data
ccRCC <- NormalizeData(ccRCC, verbose = FALSE)

ccRCC <- FindVariableFeatures(ccRCC, selection.method = "vst", nfeatures = 2000,
                              verbose = FALSE)

ccRCC <- ScaleData(ccRCC, features = row.names(ccRCC))

ccRCC <- RunPCA(ccRCC, features = VariableFeatures(object = ccRCC))

# Prepare lists for consolidating cell types
list_CD8T <- c("CD8A+ Tissue-resident", "CD8A+ NK-like",
               "CD8A+ Exhausted", "CD8A+ Exhausted IEG",
               "CD8A+ Proliferating")

list_CD4T <- c("CD4+ Activated IEG", "CD4+ Effector", 
               "CD4+ Treg", "CD4+ Naive", 
               "CD4+ Proliferating")

list_monocyte <- c("CD14+/CD16+ Monocyte", "CD14+ Monocyte")

list_NK <- c("NK HSP+", "Conventional NK")

list_macrophage <- c("TAM HLAhi", "TAM HLAint",
                     "TAM ISGint", "TAM/TCR (Ambiguos)", "TAM ISGhi")

list_DC <- c("cDC2", "pDC", 
             "cDC1")

# Consolidate cell types
ccRCC_meta <- ccRCC@meta.data

ccRCC_meta$cell_type_consolidate <- ccRCC_meta$cell_subtype

ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_CD8T] <- "CD8+ T-cell"
ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_CD4T] <- "CD4+ T-cell"
ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_monocyte] <- "Monocyte"
ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_NK] <- "NK cell"
ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_macrophage] <- "Macrophage"
ccRCC_meta$cell_type_consolidate[ccRCC_meta$cell_type_consolidate %in% list_DC] <- "DC"

# Replace old metadata with new metadata with consolidated cell types
ccRCC@meta.data <- ccRCC_meta

# Save preprocess ccRCC file
saveRDS(ccRCC, paste0(base_output_path, "ccRCC_tissue_and_PBMC_preprocessed.rds"))

####### Extract ground truth from label transfer results

##### 1) ccRCC without paired PBMC
# Get ccRCC without_paied_PBMC
ccRCC_without_paired_PBMC <- subset(ccRCC, subset = type != "PBMC")

# Check if all samples in label transfer result is in ccRCC without paired PBMC
sum(row.names(CD4T_predictions) %in% colnames(ccRCC_without_paired_PBMC)) ==
nrow(CD4T_predictions)

sum(row.names(CD8T_predictions) %in% colnames(ccRCC_without_paired_PBMC)) ==
nrow(CD8T_predictions)

sum(row.names(monocyte_predictions) %in% colnames(ccRCC_without_paired_PBMC)) ==
nrow(monocyte_predictions)

sum(row.names(NK_predictions) %in% colnames(ccRCC_without_paired_PBMC)) ==
nrow(NK_predictions)

sum(row.names(Bcell_predictions) %in% colnames(ccRCC_without_paired_PBMC)) ==
nrow(Bcell_predictions)

# Get sample IDs of potential blood immune cells
## From CD4T
CD4T_predictions_blood <- subset(CD4T_predictions, subset = predicted.id == "PBMC_CD4+ T-cell")
list_samples_blood_CD4T <- row.names(CD4T_predictions_blood)

## From CD8T
CD8T_predictions_blood <- subset(CD8T_predictions, subset = predicted.id == "PBMC_CD8+ T-cell")
list_samples_blood_CD8T <- row.names(CD8T_predictions_blood)

## From monocyte
monocyte_predictions_blood <- subset(monocyte_predictions, subset = predicted.id == "PBMC_Monocyte")
list_samples_monocyte <- row.names(monocyte_predictions_blood)

## From NK cell
NK_predictions_blood <- subset(NK_predictions, subset = predicted.id == "PBMC_NK cell")
list_samples_NK <- row.names(NK_predictions_blood)

## From B cell
Bcell_predictions_blood <- subset(Bcell_predictions, subset = predicted.id == "PBMC_B cell")
list_samples_Bcell <- row.names(Bcell_predictions_blood)

# Assign blood and tissue-resident immune cell groundtruth based on label transfer results
ccRCC_without_paired_PBMC_meta <- ccRCC_without_paired_PBMC@meta.data

ccRCC_without_paired_PBMC_meta$TB_groundtruth <- "Tissue"

ccRCC_without_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_without_paired_PBMC_meta) %in% list_samples_blood_CD4T] <- "Blood"
ccRCC_without_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_without_paired_PBMC_meta) %in% list_samples_blood_CD8T] <- "Blood"
ccRCC_without_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_without_paired_PBMC_meta) %in% list_samples_monocyte] <- "Blood"
ccRCC_without_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_without_paired_PBMC_meta) %in% list_samples_NK] <- "Blood"
ccRCC_without_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_without_paired_PBMC_meta) %in% list_samples_Bcell] <- "Blood"

# Replace old metadata with a new metadata with groudntruth
ccRCC_without_paired_PBMC@meta.data <- ccRCC_without_paired_PBMC_meta

# Save RDS file
saveRDS(ccRCC_without_paired_PBMC, 
        paste0(base_output_path ,"krishna_2021_complete_response_without_paired_PBMC_with_groudtruth.rds"))

##### 2) ccRCC with paired PBMC
# Get ccRCC with  paired PBMC
ccRCC_with_paired_PBMC <- ccRCC

# Assign blood and tissue-resident immune cell groundtruth based on label transfer results
ccRCC_with_paired_PBMC_meta <- ccRCC_with_paired_PBMC@meta.data

ccRCC_with_paired_PBMC_meta$TB_groundtruth <- "Tissue"

ccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_with_paired_PBMC_meta) %in% list_samples_blood_CD4T] <- "Blood"
ccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_with_paired_PBMC_meta) %in% list_samples_blood_CD8T] <- "Blood"
ccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_with_paired_PBMC_meta) %in% list_samples_monocyte] <- "Blood"
ccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_with_paired_PBMC_meta) %in% list_samples_NK] <- "Blood"
ccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(ccRCC_with_paired_PBMC_meta) %in% list_samples_Bcell] <- "Blood"

ccRCC_with_paired_PBMC_meta$TB_groundtruth[ccRCC_with_paired_PBMC_meta$type == "PBMC"] <- "Blood"

# Replace old metadata with a new metadata with groudntruth
ccRCC_with_paired_PBMC@meta.data <- ccRCC_with_paired_PBMC_meta

# Save RDS file
saveRDS(ccRCC_with_paired_PBMC, 
        paste0(base_output_path ,"krishna_2021_complete_response_with_paired_PBMC_with_groudtruth.rds"))

##### 3) cleaned ccRCC with paired PBMC
# Remove inferred blood contamination from tissue
cleanedccRCC_with_paired_PBMC <- ccRCC

# Assign blood and tissue-resident immune cell groundtruth based on label transfer results
cleanedccRCC_with_paired_PBMC_meta <- cleanedccRCC_with_paired_PBMC@meta.data

cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth <- "Tissue"

cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(cleanedccRCC_with_paired_PBMC_meta) %in% list_samples_blood_CD4T] <- "Blood contamination"
cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(cleanedccRCC_with_paired_PBMC_meta) %in% list_samples_blood_CD8T] <- "Blood contamination"
cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(cleanedccRCC_with_paired_PBMC_meta) %in% list_samples_monocyte] <- "Blood contamination"
cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(cleanedccRCC_with_paired_PBMC_meta) %in% list_samples_NK] <- "Blood contamination"
cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[row.names(cleanedccRCC_with_paired_PBMC_meta) %in% list_samples_Bcell] <- "Blood contamination"

cleanedccRCC_with_paired_PBMC_meta$TB_groundtruth[cleanedccRCC_with_paired_PBMC_meta$type == "PBMC"] <- "Blood"

# Replace old metadata with a new metadata with groudntruth
cleanedccRCC_with_paired_PBMC@meta.data <- cleanedccRCC_with_paired_PBMC_meta

# Remove blood contamination
cleanedccRCC_with_paired_PBMC <- subset(cleanedccRCC_with_paired_PBMC,
                                        subset = TB_groundtruth != "Blood contamination")

# Save RDS file
saveRDS(cleanedccRCC_with_paired_PBMC, 
        paste0(base_output_path ,"krishna_2021_complete_response_with_paired_PBMC_removeBloodContamination_with_groudtruth.rds"))


