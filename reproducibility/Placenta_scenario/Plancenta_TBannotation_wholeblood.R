# This file is used to perform TB annotation in placenta using the whole blood reference
# Plus, annotate subtypes of identified tissue and blood-resident immune cells

set.seed(20211122)
library(ggplot2);
library(data.table);
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library("RColorBrewer")
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

library(EnsDb.Hsapiens.v79)

setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

library(clustree)

# library(usethis)
# usethis::edit_r_environ()

'%!in%' <- function(x,y)!('%in%'(x,y))

# Read data
load("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/PL.integrated_SNG_annotated_4_27_2022.RData")


data <- PL.integrated_3
rm(PL.integrated_3)

list_hongkongs_data_name <- c("PE1", "PE2", "PE3", "PE4", "PN1", "PN2", "PN3C", "PN4C")

data.hongkong <- subset(data, subset = orig.ident %in% list_hongkongs_data_name)

# Filter to get immune cells
data.hongkong.immune <- subset(data.hongkong, subset = compartment == "Immune")

# Get the whole blood scRNA-seq reference data from Tabular Sapiens
# NOTE: Use the version with the modified gene names and got the duplicated gene IDs removed
tabular_sapiens_blood_editedGeneNames_noduplicatedGenes <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/Tabular_sapiens_blood_editedGeneNames_noduplicatedGenes_withLevel1Annotation.rds")

####### Perform blood-tissue annotation on T-cells, monocyte, NK cell
data.forQuery.immune <- data.hongkong.immune

######## NOTE: This step is run on Server due to high computational cost
data_list <- c(tabular_sapiens_blood_editedGeneNames_noduplicatedGenes, data_forQuery.immune)

wholeblood_data.anchor <- FindIntegrationAnchors(object.list = data_list, dims = 1:30, anchor.features = 3000)
wholeblood.integrated <- IntegrateData(anchorset = wholeblood_data.anchor, dims = 1:30)

wholeblood.integrated <- ScaleData(wholeblood.integrated, verbose = FALSE)
wholeblood.integrated <- RunPCA(wholeblood.integrated, npcs = 30, verbose = FALSE)
wholeblood.integrated <- RunUMAP(wholeblood.integrated, reduction = "pca", dims = 1:30, n.components = 10)
wholeblood.integrated <- FindNeighbors(wholeblood.integrated, reduction = "pca", dims = 1:30)
wholeblood.integrated <- FindClusters(wholeblood.integrated, resolution = 0.4)

# saveRDS(wholeblood.integrated, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/placenta_wholeblood_immune_integrated_forTBannotation.rds")

# Perform tissue and blood immune cells identification
wholeblood.integrated <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/placenta_wholeblood_immune_integrated_forTBannotation.rds")

# Define query and reference on common immune cell types
# Placenta side: T-cells, Monocyte, Natural Killer cells, Macrophage (HB)
## Reference side: T-cell, monocyte, nk cell, macrophage

list_common_immune_query <- c("T-cells", "Monocyte", "Natural Killer cells", "Macrophage (HB)")
list_common_immune_reference <- c("T-cell", "monocyte", "nk cell", "macrophage")

wholeblood.integrated_meta <- wholeblood.integrated@meta.data

wholeblood.integrated_meta$annotation_query_ref <- rep(NA, nrow(wholeblood.integrated_meta))

wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "T-cells"] <- "Query T-cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Monocyte"] <- "Query monocyte"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Natural Killer cells"] <- "Query NK cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Macrophage (HB)"] <- "Query macrophage"

wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "T-cell"] <- "Ref T-cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "monocyte"] <- "Ref monocyte"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "nk cell"] <- "Ref NK cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "macrophage"] <- "Ref macrophage"


wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "T-cells"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Monocyte"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Natural Killer cells"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Macrophage (HB)"] <- "Query"

wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "T-cell"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "monocyte"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "nk cell"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "macrophage"] <- "Ref"


wholeblood.integrated.common.immune <- wholeblood.integrated

wholeblood.integrated.common.immune@meta.data <-  wholeblood.integrated_meta

list_annotation_query_ref <- unlist(unique(wholeblood.integrated.common.immune$annotation_query_ref))
list_annotation_query_ref <- list_annotation_query_ref[is.na(list_annotation_query_ref) == FALSE]

wholeblood.integrated.common.immune <- subset(wholeblood.integrated.common.immune,
                                              subset = annotation_query_ref %in% list_annotation_query_ref)

# Extract umap embedding 1-10
wholeblood.integrated.common.immune_subset_xy <- wholeblood.integrated.common.immune[["umap"]]@cell.embeddings
wholeblood.integrated.common.immune_subset_meta <- wholeblood.integrated.common.immune[[]]
wholeblood.integrated.common.immune_subset_meta <- wholeblood.integrated.common.immune_subset_meta[, c("orig.ident", "annotation_query_ref")]
wholeblood.integrated.common.immune_subset_meta <- wholeblood.integrated.common.immune_subset_meta[row.names(wholeblood.integrated.common.immune_subset_meta), ]
table(row.names(wholeblood.integrated.common.immune_subset_meta) == row.names(wholeblood.integrated.common.immune_subset_xy))
wholeblood.integrated.common.immune_subset_meta_2 <- cbind(wholeblood.integrated.common.immune_subset_meta, wholeblood.integrated.common.immune_subset_xy)

###plot distance matrix and kemans cluster based on umap 1-10
library(RColorBrewer)
my.break1 <- seq(0, 15)
my.color1 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(my.break1))

my.break2 <- seq(0, 20)
my.color2 <- colorRampPalette(rev(brewer.pal(n = 14, name = "RdYlBu")))(length(my.break1))

my.break3 <- seq(0, 15)
my.color3 <- colorRampPalette(rev(brewer.pal(n = 4, name = "RdYlBu")))(length(my.break1))


# Perform blood and tissue-resident immune cell separation on common cell types
list_annotation_query_ref

DimPlot(wholeblood.integrated.common.immune, group.by = "annotation_query_ref", split.by = "orig.ident")

#############
ct = "T-cell"
temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -7.1 < UMAP 1 < 2.90; -3.80 < UMAP 2 < 7.20
temp1 <- temp1[temp1$UMAP_1 > -7.1 & temp1$UMAP_1 < 2.90 & temp1$UMAP_2 > -3.80 & temp1$UMAP_2 < 7.20, ]
ref_cid <- row.names(temp1[temp1$orig.ident == "Ref", ])
query_cid <- row.names(temp1[temp1$orig.ident == "Query", ])
temp2 <- as.matrix(temp1[, 3:12])
temp3 <- as.matrix(dist(temp2))
temp4 <- temp3[query_cid, ref_cid]

km.res <- kmeans(temp4, 2, nstart = 25)
km.res2 <- cbind(temp4, cluster = km.res$cluster)
km.res3 <- km.res2[, "cluster", drop = F]
km.res3 <- as.data.frame(km.res3)
km.res3 <- km.res3[order(km.res3$cluster), , drop = F]
km.res3$cluster <- as.factor(km.res3$cluster)
temp4 <- temp4[row.names(km.res3), ]

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_Tcells_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_Tcells_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_Tcells_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Tissue"
km.res3$origin_tb[km.res3$cluster == 2] <- "Blood"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_Tcells_UMAP_d10_query_ref_dist.RData")

#############
ct = "monocyte"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -6.4 < UMAP 1 < 2.94; -4.00 < UMAP 2 < 6.65
temp1 <- temp1[temp1$UMAP_1 > -6.4 & temp1$UMAP_1 < 2.94 & temp1$UMAP_2 > -4.00 & temp1$UMAP_2 < 6.65, ]
ref_cid <- row.names(temp1[temp1$orig.ident == "Ref", ])
query_cid <- row.names(temp1[temp1$orig.ident == "Query", ])
temp2 <- as.matrix(temp1[, 3:12])
temp3 <- as.matrix(dist(temp2))
temp4 <- temp3[query_cid, ref_cid]

km.res <- kmeans(temp4, 2, nstart = 25)
km.res2 <- cbind(temp4, cluster = km.res$cluster)
km.res3 <- km.res2[, "cluster", drop = F]
km.res3 <- as.data.frame(km.res3)
km.res3 <- km.res3[order(km.res3$cluster), , drop = F]
km.res3$cluster <- as.factor(km.res3$cluster)
temp4 <- temp4[row.names(km.res3), ]

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_monocyte_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_monocyte_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_monocyte_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Tissue"
km.res3$origin_tb[km.res3$cluster == 2] <- "Blood"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_monocyte_UMAP_d10_query_ref_dist.RData")

#############
ct = "NK cell"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: 0.07 < UMAP 1 < 2.80; -2.90 < UMAP 2 < 6.76
temp1 <- temp1[temp1$UMAP_1 > 0.07 & temp1$UMAP_1 < 2.80 & temp1$UMAP_2 > -2.90 & temp1$UMAP_2 < 6.76, ]
ref_cid <- row.names(temp1[temp1$orig.ident == "Ref", ])
query_cid <- row.names(temp1[temp1$orig.ident == "Query", ])
temp2 <- as.matrix(temp1[, 3:12])
temp3 <- as.matrix(dist(temp2))
temp4 <- temp3[query_cid, ref_cid]

km.res <- kmeans(temp4, 2, nstart = 25)
km.res2 <- cbind(temp4, cluster = km.res$cluster)
km.res3 <- km.res2[, "cluster", drop = F]
km.res3 <- as.data.frame(km.res3)
km.res3 <- km.res3[order(km.res3$cluster), , drop = F]
km.res3$cluster <- as.factor(km.res3$cluster)
temp4 <- temp4[row.names(km.res3), ]

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_NKcell_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_NKcell_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_NKcell_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_NKcell_UMAP_d10_query_ref_dist.RData")

#############
ct = "macrophage"
temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -3.93 < UMAP 1 < 2.65; -3.73 < UMAP 2 < 5.66
temp1 <- temp1[temp1$UMAP_1 > -3.93 & temp1$UMAP_1 < 2.65 & temp1$UMAP_2 > -3.73 & temp1$UMAP_2 < 5.66, ]
ref_cid <- row.names(temp1[temp1$orig.ident == "Ref", ])
query_cid <- row.names(temp1[temp1$orig.ident == "Query", ])
temp2 <- as.matrix(temp1[, 3:12])
temp3 <- as.matrix(dist(temp2))
temp4 <- temp3[query_cid, ref_cid]

km.res <- kmeans(temp4, 2, nstart = 25)
km.res2 <- cbind(temp4, cluster = km.res$cluster)
km.res3 <- km.res2[, "cluster", drop = F]
km.res3 <- as.data.frame(km.res3)
km.res3 <- km.res3[order(km.res3$cluster), , drop = F]
km.res3$cluster <- as.factor(km.res3$cluster)
temp4 <- temp4[row.names(km.res3), ]

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_macrophage_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_macrophage_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_macrophage_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_macrophage_UMAP_d10_query_ref_dist.RData")

##collect all data
temp_meta <- NULL
load("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_Tcells_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_monocyte_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_NKcell_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Wholeblood_macrophage_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)

save(temp_meta, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/tissue_blood_annotation_usingWholeBlood.RData")

### pull results and annotate tissue and blood origin
data_TB_annotation <- data.hongkong

data_TB_annotation_meta <- data_TB_annotation@meta.data

## tissue and blood origin assignment
data_TB_annotation_meta_1 <- data_TB_annotation_meta[!(row.names(data_TB_annotation_meta) %in% row.names(temp_meta)), ]

data_TB_annotation_meta_1$TB_origin_wholeblood <- rep("Undetermined", dim(data_TB_annotation_meta_1)[1])

data_TB_annotation_meta_1$TB_origin_wholeblood[data_TB_annotation_meta_1$compartment == "Trophoblast"|
                                                 data_TB_annotation_meta_1$compartment == "Stromal"] <- "Tissue"

data_TB_annotation_meta_1$TB_origin_wholeblood[data_TB_annotation_meta_1$annotation == "Erythrocyte"] <- "Blood"

data_TB_annotation_meta_2 <- data_TB_annotation_meta_1[, "TB_origin_wholeblood", drop = F]

temp_meta <- temp_meta[, "origin_tb", drop = F]
colnames(temp_meta) <- "TB_origin_wholeblood"

temp_meta <- rbind(temp_meta, data_TB_annotation_meta_2)

temp_meta <- temp_meta[row.names(data_TB_annotation_meta), , drop = F]

data_TB_annotation[["TB_origin_wholeblood"]] <- temp_meta$TB_origin_wholeblood
                    
# Save TB annotated file
# saveRDS(data_TB_annotation, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")

# Assign maternal-fetal origin and blood & tissue-resident cells
data_MF_TB_annotation <- data_TB_annotation

data_MF_TB_annotation_meta <- data_MF_TB_annotation@meta.data

data_MF_TB_annotation_meta$MF_TB_origin_refWholeblood <- rep(NA, nrow(data_MF_TB_annotation_meta))

data_MF_TB_annotation_meta$MF_TB_origin_refWholeblood[(data_MF_TB_annotation_meta$MF_Origin == "Fetal") &
                                                        (data_MF_TB_annotation_meta$TB_origin_wholeblood == "Tissue")] <- "Fetal Tissue"

data_MF_TB_annotation_meta$MF_TB_origin_refWholeblood[(data_MF_TB_annotation_meta$MF_Origin == "Fetal") &
                                                        (data_MF_TB_annotation_meta$TB_origin_wholeblood == "Blood")] <- "Fetal Blood"

data_MF_TB_annotation_meta$MF_TB_origin_refWholeblood[(data_MF_TB_annotation_meta$MF_Origin == "Maternal") &
                                                        (data_MF_TB_annotation_meta$TB_origin_wholeblood == "Tissue")] <- "Maternal Tissue"

data_MF_TB_annotation_meta$MF_TB_origin_refWholeblood[(data_MF_TB_annotation_meta$MF_Origin == "Maternal") &
                                                        (data_MF_TB_annotation_meta$TB_origin_wholeblood == "Blood")] <- "Maternal Blood"
                             
data_MF_TB_annotation@meta.data <- data_MF_TB_annotation_meta

saveRDS(data_MF_TB_annotation, file = "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")

##############  Visualization

data_MF_TB_annotation_MB <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Blood")
data_MF_TB_annotation_MT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Tissue")
data_MF_TB_annotation_FB <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Blood")
data_MF_TB_annotation_FT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Tissue")


#### UMAP of all cell gropued by origins

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/UMAP_MF_TB_refWholeBlood.pdf", width = 17, height = 12)

DimPlot(data_MF_TB_annotation, group.by = "annotation", split.by = "MF_TB_origin_refWholeblood", ncol = 2)
DimPlot(data_MF_TB_annotation, group.by = "annotation", split.by = "MF_TB_origin_refWholeblood", ncol = 2, label = TRUE)

DimPlot(data_MF_TB_annotation_MB, group.by = "annotation")
DimPlot(data_MF_TB_annotation_MB, group.by = "annotation", label = TRUE)

DimPlot(data_MF_TB_annotation_MT, group.by = "annotation")
DimPlot(data_MF_TB_annotation_MT, group.by = "annotation", label = TRUE)

DimPlot(data_MF_TB_annotation_FB, group.by = "annotation")
DimPlot(data_MF_TB_annotation_FB, group.by = "annotation", label = TRUE)

DimPlot(data_MF_TB_annotation_FT, group.by = "annotation")
DimPlot(data_MF_TB_annotation_FT, group.by = "annotation", label = TRUE)

dev.off()

save.image("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TBannotation_wholeblood.RData")

