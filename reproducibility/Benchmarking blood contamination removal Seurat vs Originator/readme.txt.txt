###############################################################
#### FILE: derive_groundtruth_from_tissue_Krishna_et_al_2021.R
#### DESCRIPTION:
# This file is used to derive ground truth from PBMC and impure tissue data using Seurat label transfer
# Dataset in this study: Krishna et al. (2021) - Kidney
#   - Use a patient with complete response as an example
#   - Path to file: /nfs/dcmb-lgarmire/yangiwen/workspace/originator/data/krishna/3ca/complete_response_all.rds
# Given that each patient has paired PBMC and ccRCC tissues
#   1) Use label transfer to identify PBMC in ccRCC tissues by using
#       PBMC as "reference" and ccRCC tissue as "query"
#   2) Remove identified blood immune cells in ccRCC tissues as and
#       treat the rest immune cells in CLEANED ccRCC as "tissue-resident immume cells" with paired PBMC as "blood immune cells"
#   3) Use these labels as GROUND TRUTH for Originator
#
# NOTE: Use Seurat version 5.0.1

###############################################################
#### FILE: assign_groundtruth_from_tissue_Krishna_et_al_2021.R
#### DESCRIPTION:
# This file is used to assign ground truth derived from derive_groundtruth_from_tissue_Krishna_et_al_2021.R to prepare dataset for
#	running Originator

###############################################################
#### FILE: assign_groundtruth_from_tissue_Krishna_et_al_2021.R
#### DESCRIPTION:
# This file is used to assign groundtruth derived from Originator using wholeblood as reference and paired ccRCC as query
# Blood contamination in ccRCC will be removed. The rest will be considered as tissue-resident immune cells. Paired PBMC will be 
#	considered as ground truth for blood immune cells

###############################################################
#### FILE: visualizeF1_TB_annotation_krishna_et_al_2021_seurat_vs_originator.R
#### DESCRIPTION:
# This file is used to plot a box plot of F1 score comparing Originator-based and
#	Seurat-based in  TB annotation of cleaned ccRCC & paired PBMC



