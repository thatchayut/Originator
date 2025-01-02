# This file is used to assign groundtruth derived from Originator using
#   wholeblood as reference and paired ccRCC as query
# Blood contamination in ccRCC will be removed. The rest will be considered as tissue-resident immune cells. Paired PBMC will be considered as ground truth
#   for blood immune cells

##### NOTE: Use Seurat version 4.3.0
library(Seurat)

setwd("~/temp/soft_link_nfs_originator_GB_revision")

# Specify path to output folder
base_output_path <- "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/"

# Read input file
ccRCC_noPBMC <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/krishna_2021_complete_response_without_paired_PBMC.rds")

########### OPTION 1: Use paired PBMC as a reference for Originator to identify contamination

# Read prediction files performed by originator
CD4T_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator/cd4/0/CD4_annotated.rds")

CD8T_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator/cd8/0/CD8_annotated.rds")

monocyte_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator/mono/0/Monocyte_annotated.rds")

NK_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator/nk/0/NK-cell_annotated.rds")

Bcell_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator/b/0/B-cell_annotated.rds")

# Get list of samples considered as blood contamination
## From CD4T
CD4T_predictions_blood <- subset(CD4T_predictions, subset = origin_tb == "Blood")
list_samples_blood_CD4T <- row.names(CD4T_predictions_blood@meta.data)

## From CD8T
CD8T_predictions_blood <- subset(CD8T_predictions, subset = origin_tb == "Blood")
list_samples_blood_CD8T <- row.names(CD8T_predictions_blood@meta.data)

## From monocyte
monocyte_predictions_blood <- subset(monocyte_predictions, subset = origin_tb == "Blood")
list_samples_monocyte <- row.names(monocyte_predictions_blood@meta.data)

## From NK cell
NK_predictions_blood <- subset(NK_predictions, subset = origin_tb == "Blood")
list_samples_NK <- row.names(NK_predictions_blood@meta.data)

## From B cell
Bcell_predictions_blood <- subset(Bcell_predictions, subset = origin_tb == "Blood")
list_samples_Bcell <- row.names(Bcell_predictions_blood@meta.data)

# Remove inferred blood contamination from tissue
ccRCC_with_paired_PBMC <- readRDS("/nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds")

cleanedccRCC_with_paired_PBMC <- ccRCC_with_paired_PBMC

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
        paste0(base_output_path ,"krishna_2021_complete_response_with_paired_PBMC_removeBloodContaminationByOriginator_with_groudtruth.rds"))

########### OPTION 2: Use wholeblood as a reference for Originator to identify contamination
# Read prediction files performed by originator
CD4T_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator_wholebloodRef/cd4/0/CD4_annotated.rds")

CD8T_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator_wholebloodRef/cd8/0/CD8_annotated.rds")

monocyte_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator_wholebloodRef/mono/0/Monocyte_annotated.rds")

NK_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator_wholebloodRef/nk/0/NK-cell_annotated.rds")

Bcell_predictions <- readRDS("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/results/kidney_krishna_seeds/deriveGroundtruthByOriginator_wholebloodRef/b/0/B-cell_annotated.rds")

# Get list of samples considered as blood contamination
## From CD4T
CD4T_predictions_blood <- subset(CD4T_predictions, subset = origin_tb == "Blood")
list_samples_blood_CD4T <- row.names(CD4T_predictions_blood@meta.data)

## From CD8T
CD8T_predictions_blood <- subset(CD8T_predictions, subset = origin_tb == "Blood")
list_samples_blood_CD8T <- row.names(CD8T_predictions_blood@meta.data)

## From monocyte
monocyte_predictions_blood <- subset(monocyte_predictions, subset = origin_tb == "Blood")
list_samples_monocyte <- row.names(monocyte_predictions_blood@meta.data)

## From NK cell
NK_predictions_blood <- subset(NK_predictions, subset = origin_tb == "Blood")
list_samples_NK <- row.names(NK_predictions_blood@meta.data)

## From B cell
Bcell_predictions_blood <- subset(Bcell_predictions, subset = origin_tb == "Blood")
list_samples_Bcell <- row.names(Bcell_predictions_blood@meta.data)

# Remove inferred blood contamination from tissue
ccRCC_with_paired_PBMC <- readRDS("/nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds")

cleanedccRCC_with_paired_PBMC <- ccRCC_with_paired_PBMC

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
        paste0(base_output_path ,"krishna_2021_complete_response_with_paired_PBMC_removeBloodContaminationByOriginator_wholebloodRef_with_groudtruth.rds"))

