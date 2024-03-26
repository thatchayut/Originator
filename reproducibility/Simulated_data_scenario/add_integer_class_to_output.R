# This file is used to add integer class to the dataframe output to be used when evaluating the performance


library(dplyr)

# Read input data
## iteration 1
df_result_iter1_tcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell.csv")
df_result_iter1_monocyte <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte.csv")
df_result_iter1_bcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell.csv")
df_result_iter1_overall <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune.csv")

## iteration 2
df_result_iter2_tcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell.csv")
df_result_iter2_monocyte <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte.csv")
df_result_iter2_bcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell.csv")
df_result_iter2_overall <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune.csv")

## iteration 3
df_result_iter3_tcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell.csv")
df_result_iter3_monocyte <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte.csv")
df_result_iter3_bcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell.csv")
df_result_iter3_overall <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune.csv")

## iteration 4
df_result_iter4_tcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell.csv")
df_result_iter4_monocyte <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte.csv")
df_result_iter4_bcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell.csv")
df_result_iter4_overall <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune.csv")

## iteration 5
df_result_iter5_tcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell.csv")
df_result_iter5_monocyte <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte.csv")
df_result_iter5_bcell <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell.csv")
df_result_iter5_overall <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune.csv")

# Add integer classes to ground truths and predictions (Tissue = 1, Blood = 0)
## iteration 1
### T cell
df_result_iter1_tcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter1_tcell))
df_result_iter1_tcell$origin_groundtruth_integer[df_result_iter1_tcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter1_tcell$origin_groundtruth_integer[df_result_iter1_tcell$origin_groundtruth == "Blood"] <- 0

df_result_iter1_tcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter1_tcell))
df_result_iter1_tcell$TB_origin_wholeblood_integer[df_result_iter1_tcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter1_tcell$TB_origin_wholeblood_integer[df_result_iter1_tcell$TB_origin_wholeblood == "Blood"] <- 0

### Monocyte
df_result_iter1_monocyte$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter1_monocyte))
df_result_iter1_monocyte$origin_groundtruth_integer[df_result_iter1_monocyte$origin_groundtruth == "Tissue"] <- 1
df_result_iter1_monocyte$origin_groundtruth_integer[df_result_iter1_monocyte$origin_groundtruth == "Blood"] <- 0

df_result_iter1_monocyte$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter1_monocyte))
df_result_iter1_monocyte$TB_origin_wholeblood_integer[df_result_iter1_monocyte$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter1_monocyte$TB_origin_wholeblood_integer[df_result_iter1_monocyte$TB_origin_wholeblood == "Blood"] <- 0

### B cell
df_result_iter1_bcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter1_bcell))
df_result_iter1_bcell$origin_groundtruth_integer[df_result_iter1_bcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter1_bcell$origin_groundtruth_integer[df_result_iter1_bcell$origin_groundtruth == "Blood"] <- 0

df_result_iter1_bcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter1_bcell))
df_result_iter1_bcell$TB_origin_wholeblood_integer[df_result_iter1_bcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter1_bcell$TB_origin_wholeblood_integer[df_result_iter1_bcell$TB_origin_wholeblood == "Blood"] <- 0

### Overall 
df_result_iter1_overall$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter1_overall))
df_result_iter1_overall$origin_groundtruth_integer[df_result_iter1_overall$origin_groundtruth == "Tissue"] <- 1
df_result_iter1_overall$origin_groundtruth_integer[df_result_iter1_overall$origin_groundtruth == "Blood"] <- 0

df_result_iter1_overall$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter1_overall))
df_result_iter1_overall$TB_origin_wholeblood_integer[df_result_iter1_overall$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter1_overall$TB_origin_wholeblood_integer[df_result_iter1_overall$TB_origin_wholeblood == "Blood"] <- 0

## iteration 2
### T cell
df_result_iter2_tcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter2_tcell))
df_result_iter2_tcell$origin_groundtruth_integer[df_result_iter2_tcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter2_tcell$origin_groundtruth_integer[df_result_iter2_tcell$origin_groundtruth == "Blood"] <- 0

df_result_iter2_tcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter2_tcell))
df_result_iter2_tcell$TB_origin_wholeblood_integer[df_result_iter2_tcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter2_tcell$TB_origin_wholeblood_integer[df_result_iter2_tcell$TB_origin_wholeblood == "Blood"] <- 0

### Monocyte
df_result_iter2_monocyte$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter2_monocyte))
df_result_iter2_monocyte$origin_groundtruth_integer[df_result_iter2_monocyte$origin_groundtruth == "Tissue"] <- 1
df_result_iter2_monocyte$origin_groundtruth_integer[df_result_iter2_monocyte$origin_groundtruth == "Blood"] <- 0

df_result_iter2_monocyte$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter2_monocyte))
df_result_iter2_monocyte$TB_origin_wholeblood_integer[df_result_iter2_monocyte$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter2_monocyte$TB_origin_wholeblood_integer[df_result_iter2_monocyte$TB_origin_wholeblood == "Blood"] <- 0

### B cell
df_result_iter2_bcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter2_bcell))
df_result_iter2_bcell$origin_groundtruth_integer[df_result_iter2_bcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter2_bcell$origin_groundtruth_integer[df_result_iter2_bcell$origin_groundtruth == "Blood"] <- 0

df_result_iter2_bcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter2_bcell))
df_result_iter2_bcell$TB_origin_wholeblood_integer[df_result_iter2_bcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter2_bcell$TB_origin_wholeblood_integer[df_result_iter2_bcell$TB_origin_wholeblood == "Blood"] <- 0

### Overall 
df_result_iter2_overall$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter2_overall))
df_result_iter2_overall$origin_groundtruth_integer[df_result_iter2_overall$origin_groundtruth == "Tissue"] <- 1
df_result_iter2_overall$origin_groundtruth_integer[df_result_iter2_overall$origin_groundtruth == "Blood"] <- 0

df_result_iter2_overall$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter2_overall))
df_result_iter2_overall$TB_origin_wholeblood_integer[df_result_iter2_overall$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter2_overall$TB_origin_wholeblood_integer[df_result_iter2_overall$TB_origin_wholeblood == "Blood"] <- 0

## iteration 3
### T cell
df_result_iter3_tcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter3_tcell))
df_result_iter3_tcell$origin_groundtruth_integer[df_result_iter3_tcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter3_tcell$origin_groundtruth_integer[df_result_iter3_tcell$origin_groundtruth == "Blood"] <- 0

df_result_iter3_tcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter3_tcell))
df_result_iter3_tcell$TB_origin_wholeblood_integer[df_result_iter3_tcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter3_tcell$TB_origin_wholeblood_integer[df_result_iter3_tcell$TB_origin_wholeblood == "Blood"] <- 0

### Monocyte
df_result_iter3_monocyte$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter3_monocyte))
df_result_iter3_monocyte$origin_groundtruth_integer[df_result_iter3_monocyte$origin_groundtruth == "Tissue"] <- 1
df_result_iter3_monocyte$origin_groundtruth_integer[df_result_iter3_monocyte$origin_groundtruth == "Blood"] <- 0

df_result_iter3_monocyte$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter3_monocyte))
df_result_iter3_monocyte$TB_origin_wholeblood_integer[df_result_iter3_monocyte$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter3_monocyte$TB_origin_wholeblood_integer[df_result_iter3_monocyte$TB_origin_wholeblood == "Blood"] <- 0

### B cell
df_result_iter3_bcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter3_bcell))
df_result_iter3_bcell$origin_groundtruth_integer[df_result_iter3_bcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter3_bcell$origin_groundtruth_integer[df_result_iter3_bcell$origin_groundtruth == "Blood"] <- 0

df_result_iter3_bcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter3_bcell))
df_result_iter3_bcell$TB_origin_wholeblood_integer[df_result_iter3_bcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter3_bcell$TB_origin_wholeblood_integer[df_result_iter3_bcell$TB_origin_wholeblood == "Blood"] <- 0

### Overall 
df_result_iter3_overall$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter3_overall))
df_result_iter3_overall$origin_groundtruth_integer[df_result_iter3_overall$origin_groundtruth == "Tissue"] <- 1
df_result_iter3_overall$origin_groundtruth_integer[df_result_iter3_overall$origin_groundtruth == "Blood"] <- 0

df_result_iter3_overall$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter3_overall))
df_result_iter3_overall$TB_origin_wholeblood_integer[df_result_iter3_overall$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter3_overall$TB_origin_wholeblood_integer[df_result_iter3_overall$TB_origin_wholeblood == "Blood"] <- 0

## iteration 4
### T cell
df_result_iter4_tcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter4_tcell))
df_result_iter4_tcell$origin_groundtruth_integer[df_result_iter4_tcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter4_tcell$origin_groundtruth_integer[df_result_iter4_tcell$origin_groundtruth == "Blood"] <- 0

df_result_iter4_tcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter4_tcell))
df_result_iter4_tcell$TB_origin_wholeblood_integer[df_result_iter4_tcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter4_tcell$TB_origin_wholeblood_integer[df_result_iter4_tcell$TB_origin_wholeblood == "Blood"] <- 0

### Monocyte
df_result_iter4_monocyte$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter4_monocyte))
df_result_iter4_monocyte$origin_groundtruth_integer[df_result_iter4_monocyte$origin_groundtruth == "Tissue"] <- 1
df_result_iter4_monocyte$origin_groundtruth_integer[df_result_iter4_monocyte$origin_groundtruth == "Blood"] <- 0

df_result_iter4_monocyte$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter4_monocyte))
df_result_iter4_monocyte$TB_origin_wholeblood_integer[df_result_iter4_monocyte$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter4_monocyte$TB_origin_wholeblood_integer[df_result_iter4_monocyte$TB_origin_wholeblood == "Blood"] <- 0

### B cell
df_result_iter4_bcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter4_bcell))
df_result_iter4_bcell$origin_groundtruth_integer[df_result_iter4_bcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter4_bcell$origin_groundtruth_integer[df_result_iter4_bcell$origin_groundtruth == "Blood"] <- 0

df_result_iter4_bcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter4_bcell))
df_result_iter4_bcell$TB_origin_wholeblood_integer[df_result_iter4_bcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter4_bcell$TB_origin_wholeblood_integer[df_result_iter4_bcell$TB_origin_wholeblood == "Blood"] <- 0

### Overall 
df_result_iter4_overall$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter4_overall))
df_result_iter4_overall$origin_groundtruth_integer[df_result_iter4_overall$origin_groundtruth == "Tissue"] <- 1
df_result_iter4_overall$origin_groundtruth_integer[df_result_iter4_overall$origin_groundtruth == "Blood"] <- 0

df_result_iter4_overall$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter4_overall))
df_result_iter4_overall$TB_origin_wholeblood_integer[df_result_iter4_overall$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter4_overall$TB_origin_wholeblood_integer[df_result_iter4_overall$TB_origin_wholeblood == "Blood"] <- 0

## iteration 5
### T cell
df_result_iter5_tcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter5_tcell))
df_result_iter5_tcell$origin_groundtruth_integer[df_result_iter5_tcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter5_tcell$origin_groundtruth_integer[df_result_iter5_tcell$origin_groundtruth == "Blood"] <- 0

df_result_iter5_tcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter5_tcell))
df_result_iter5_tcell$TB_origin_wholeblood_integer[df_result_iter5_tcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter5_tcell$TB_origin_wholeblood_integer[df_result_iter5_tcell$TB_origin_wholeblood == "Blood"] <- 0

### Monocyte
df_result_iter5_monocyte$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter5_monocyte))
df_result_iter5_monocyte$origin_groundtruth_integer[df_result_iter5_monocyte$origin_groundtruth == "Tissue"] <- 1
df_result_iter5_monocyte$origin_groundtruth_integer[df_result_iter5_monocyte$origin_groundtruth == "Blood"] <- 0

df_result_iter5_monocyte$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter5_monocyte))
df_result_iter5_monocyte$TB_origin_wholeblood_integer[df_result_iter5_monocyte$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter5_monocyte$TB_origin_wholeblood_integer[df_result_iter5_monocyte$TB_origin_wholeblood == "Blood"] <- 0

### B cell
df_result_iter5_bcell$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter5_bcell))
df_result_iter5_bcell$origin_groundtruth_integer[df_result_iter5_bcell$origin_groundtruth == "Tissue"] <- 1
df_result_iter5_bcell$origin_groundtruth_integer[df_result_iter5_bcell$origin_groundtruth == "Blood"] <- 0

df_result_iter5_bcell$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter5_bcell))
df_result_iter5_bcell$TB_origin_wholeblood_integer[df_result_iter5_bcell$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter5_bcell$TB_origin_wholeblood_integer[df_result_iter5_bcell$TB_origin_wholeblood == "Blood"] <- 0

### Overall 
df_result_iter5_overall$origin_groundtruth_integer <- rep(NA, nrow(df_result_iter5_overall))
df_result_iter5_overall$origin_groundtruth_integer[df_result_iter5_overall$origin_groundtruth == "Tissue"] <- 1
df_result_iter5_overall$origin_groundtruth_integer[df_result_iter5_overall$origin_groundtruth == "Blood"] <- 0

df_result_iter5_overall$TB_origin_wholeblood_integer <- rep(NA, nrow(df_result_iter5_overall))
df_result_iter5_overall$TB_origin_wholeblood_integer[df_result_iter5_overall$TB_origin_wholeblood == "Tissue"] <- 1
df_result_iter5_overall$TB_origin_wholeblood_integer[df_result_iter5_overall$TB_origin_wholeblood == "Blood"] <- 0


# write output
## iteration 1
write.csv(df_result_iter1_tcell,
           "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter1_monocyte, 
          "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter1_bcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter1_overall, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_1/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv",
          row.names = FALSE)

## iteration 2
write.csv(df_result_iter2_tcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter2_monocyte, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter2_bcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter2_overall, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_2/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv",
          row.names = FALSE)

## iteration 3
write.csv(df_result_iter3_tcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter3_monocyte, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter3_bcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter3_overall, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_3/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv",
          row.names = FALSE)

## iteration 4
write.csv(df_result_iter4_tcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter4_monocyte, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter4_bcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter4_overall, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_4/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv",
          row.names = FALSE)

## iteration 5
write.csv(df_result_iter5_tcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_tcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter5_monocyte, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_monocyte_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter5_bcell, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_Bcell_integerClasses.csv",
          row.names = FALSE)
write.csv(df_result_iter5_overall, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/artificial_tissue_blood_mixing_wholeblood_ref/results/iteration_5/meta_data_merged_tissue_blood_residents_DEMO_TB_TuN_annotation_usingWholeblood_commonImmune_integerClasses.csv",
          row.names = FALSE)


