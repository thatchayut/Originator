# This file is used to derive ground truth from PBMC and impure tissue data
# Dataset in this study: Krishna et al. (2021) - Kidney
#   - Use a patient with complete response as an example
#   - Path to file: /nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds
# Given that each patient has paired PBMC and ccRCC tissues
#   1) Use label transfer to identify PBMC in ccRCC tissues by using
#       PBMC as "reference" and ccRCC tissue as "query"
#   2) Remove identified blood immune cells in ccRCC tissues as and
#       treat the rest immune cells in CLEANED ccRCC as "tissue-resident immume cells" with paired PBMC as "blood immune cells"
#   3) Use these labels as GROUND TRUTH for Originator

##### NOTE: Use Seurat version 5.0.1

library(Seurat)

setwd("~/temp/soft_link_nfs_originator_GB_revision")

# Specify path to output folder
base_output_path <- "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/"

# Read input file
ccRCC <- readRDS("/nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds")

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

# Get ccRCC without PBMC
ccRCC_noPBMC <- subset(ccRCC, subset = type != "PBMC")

# add cell type to column "cell_type_consolidate_for_query"
ccRCC_noPBMC$cell_type_consolidate_for_query <- ccRCC_noPBMC$cell_type_consolidate

saveRDS(ccRCC_noPBMC, paste0(base_output_path, "krishna_2021_complete_response_without_paired_PBMC.rds"))

# Extract PBMC
PBMC <- subset(ccRCC, subset = type == "PBMC")

# add cell type to column "cell_type_consolidate_for_reference"
PBMC$cell_type_consolidate_for_reference <- PBMC$cell_type_consolidate

saveRDS(PBMC, paste0(base_output_path, "krishna_2021_complete_response_only_paired_PBMC.rds"))

# Extract ccRCC tissues
ccRCC_tissues <- subset(ccRCC, subset = type %in% c("Normal", "Tumor"))

# Prepare a column in PBMC data for PBMC_<cell_type> for later use as ground truth
PBMC_meta <- PBMC@meta.data
PBMC_meta$tissue_with_celltype <- paste0("PBMC_", PBMC_meta$cell_type_consolidate)

# Add metadata back to PBMC seurat object
PBMC@meta.data <- PBMC_meta

###### Perform label transfer for each PBMC cell type separately for accuracy
##### Refer to this blog: https://satijalab.org/seurat/articles/integration_mapping
#### CD4+ T-cell
query_CD4T <- subset(ccRCC_tissues, subset = cell_type_consolidate == "CD4+ T-cell")

# OPTION 1: Use only CD4+ T-cell as reference
# Extract CD4+ T-cell for query and ref
# ref_CD4T <- subset(PBMC, subset = cell_type_consolidate == "CD4+ T-cell")

# OPTION 2:Use the whole PBMC as a reference
#        only CD4+ T-cells predicted as PBMC_CD4+ T-cells will be reassigned to CD4
ref_PBMC_for_CD4T <- PBMC
ref_PBMC_for_CD4T_meta <- ref_PBMC_for_CD4T@meta.data
ref_PBMC_for_CD4T_meta$tissue_with_celltype[ref_PBMC_for_CD4T_meta$tissue_with_celltype != "PBMC_CD4+ T-cell"] <- "NOT PBMC_CD4+ T-cell"
ref_PBMC_for_CD4T@meta.data <- ref_PBMC_for_CD4T_meta

# Find integration anchors
CD4T_anchors <- FindTransferAnchors(reference = ref_PBMC_for_CD4T, query = query_CD4T,
                                    reference.reduction = "pca")

CD4T_predictions <- TransferData(anchorset = CD4T_anchors, refdata = ref_PBMC_for_CD4T$tissue_with_celltype,
                                 dims = 1:30)

#### CD8+ T-cell
query_CD8T <- subset(ccRCC_tissues, subset = cell_type_consolidate == "CD8+ T-cell")

# OPTION 1: Use only CD8+ T-cell as reference
# Extract CD8+ T-cell for query and ref
# ref_CD8T <- subset(PBMC, subset = cell_type_consolidate == "CD8+ T-cell")

# OPTION 2:Use the whole PBMC as a reference
#        only CD8+ T-cells predicted as PBMC_CD8+ T-cells will be reassigned to CD8
ref_PBMC_for_CD8T <- PBMC
ref_PBMC_for_CD8T_meta <- ref_PBMC_for_CD8T@meta.data
ref_PBMC_for_CD8T_meta$tissue_with_celltype[ref_PBMC_for_CD8T_meta$tissue_with_celltype != "PBMC_CD8+ T-cell"] <- "NOT PBMC_CD8+ T-cell"
ref_PBMC_for_CD8T@meta.data <- ref_PBMC_for_CD8T_meta

# Find integration anchors
CD8T_anchors <- FindTransferAnchors(reference = ref_PBMC_for_CD8T, query = query_CD8T,
                                    reference.reduction = "pca")

CD8T_predictions <- TransferData(anchorset = CD8T_anchors, refdata = ref_PBMC_for_CD8T$tissue_with_celltype,
                                 dims = 1:30)

#### Monocyte
query_monocyte <- subset(ccRCC_tissues, subset = cell_type_consolidate == "Monocyte")

# OPTION 1: Use only monocyte as reference
# Extract monocyte for query and ref
# ref_monocyte <- subset(PBMC, subset = cell_type_consolidate == "Monocyte")

# OPTION 2:Use the whole PBMC as a reference
#        only monocytes predicted as monocytes will be reassigned to monocytes
ref_PBMC_for_monocyte <- PBMC
ref_PBMC_for_monocyte_meta <- ref_PBMC_for_monocyte@meta.data
ref_PBMC_for_monocyte_meta$tissue_with_celltype[ref_PBMC_for_monocyte_meta$tissue_with_celltype != "PBMC_Monocyte"] <- "NOT PBMC_Monocyte"
ref_PBMC_for_monocyte@meta.data <- ref_PBMC_for_monocyte_meta

# Find integration anchors
monocyte_anchors <- FindTransferAnchors(reference = ref_PBMC_for_monocyte, query = query_monocyte,
                                    reference.reduction = "pca")

monocyte_predictions <- TransferData(anchorset = monocyte_anchors, refdata = ref_PBMC_for_monocyte$tissue_with_celltype,
                                 dims = 1:30)

#### NK cells
query_NK <- subset(ccRCC_tissues, subset = cell_type_consolidate == "NK cell")

# OPTION 1: Use only NK cell as reference
# Extract NK cell for query and ref
# ref_NK <- subset(PBMC, subset = cell_type_consolidate == "NK cell")

# OPTION 2:Use the whole PBMC as a reference
#        only NK cell predicted as NK cell will be reassigned to NK cell
ref_PBMC_for_NK <- PBMC
ref_PBMC_for_NK_meta <- ref_PBMC_for_NK@meta.data
ref_PBMC_for_NK_meta$tissue_with_celltype[ref_PBMC_for_NK_meta$tissue_with_celltype != "PBMC_NK cell"] <- "NOT PBMC_NK cell"
ref_PBMC_for_NK@meta.data <- ref_PBMC_for_NK_meta

# Find integration anchors
NK_anchors <- FindTransferAnchors(reference = ref_PBMC_for_NK, query = query_NK,
                                        reference.reduction = "pca")

NK_predictions <- TransferData(anchorset = NK_anchors, refdata = ref_PBMC_for_NK$tissue_with_celltype,
                                     dims = 1:30)

#### B cell
query_Bcell <- subset(ccRCC_tissues, subset = cell_type_consolidate == "B cell")

# OPTION 1: Use only B cell as reference
# Extract B cell for query and ref
# ref_Bcell <- subset(PBMC, subset = cell_type_consolidate == "B cell")

# OPTION 2:Use the whole PBMC as a reference
#        only B cell predicted as B cell will be reassigned to B cell
ref_PBMC_for_Bcell <- PBMC
ref_PBMC_for_Bcell_meta <- ref_PBMC_for_Bcell@meta.data
ref_PBMC_for_Bcell_meta$tissue_with_celltype[ref_PBMC_for_Bcell_meta$tissue_with_celltype != "PBMC_B cell"] <- "NOT PBMC_B cell"
ref_PBMC_for_Bcell@meta.data <- ref_PBMC_for_Bcell_meta

# Find integration anchors
Bcell_anchors <- FindTransferAnchors(reference = ref_PBMC_for_Bcell, query = query_Bcell,
                                  reference.reduction = "pca")

# k.weight is set to 27 (less than # of available anchors of 28; default = 50)
Bcell_predictions <- TransferData(anchorset = Bcell_anchors, refdata = ref_PBMC_for_Bcell$tissue_with_celltype,
                               dims = 1:30, k.weight = 27)

######## Save predictions to .rds files
saveRDS(CD4T_predictions, paste0(base_output_path, "label_transfer_result_CD4T.rds"))

saveRDS(CD8T_predictions, paste0(base_output_path, "label_transfer_result_CD8T.rds"))

saveRDS(monocyte_predictions, paste0(base_output_path, "label_transfer_result_Monocyte.rds"))

saveRDS(NK_predictions, paste0(base_output_path, "label_transfer_result_NK.rds"))

saveRDS(Bcell_predictions, paste0(base_output_path, "label_transfer_result_Bcell.rds"))









