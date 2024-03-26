# This file is used to perform TB annotation in PDAC using the whole blood reference
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

setwd("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4")

library(clustree)


# library(usethis) 
# usethis::edit_r_environ()

'%!in%' <- function(x,y)!('%in%'(x,y))


# Read data 
data <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/PDAC_TB_TuN_annotation_variant_analysis_4_addMonocytes_afterMonocyteTBannotated_020824.rds")

# Assign compartments based on the atlas that will be used as a reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10460409/
# # Compartment: Stroma, Endothelium, Myeloid, Lymphocytes
data_meta <- data@meta.data

# Perform TB annotation using the whole blood scRNA-seq reference from tabular sapiens (https://cellxgene.cziscience.com/e/983d5ec9-40e8-4512-9e65-a572a9c486cb.cxg/)
tabular_sapiens_blood <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/Tabular_sapiens_blood.rds")

# Combine substypes into level-1 annotation for TB annotation

tabular_sapiens_blood_meta <- tabular_sapiens_blood@meta.data

tabular_sapiens_blood_meta$level1_annotation <- tabular_sapiens_blood_meta$cell_ontology_class

list_tcells <- c("cd4-positive, alpha-beta memory t cell", "cd8-positive, alpha-beta cytokine secreting effector t cell",
                 "type i nk t cell", "cd8-positive, alpha-beta t cell", "t cell", 
                 "naive thymus-derived cd4-positive, alpha-beta t cell", "cd4-positive, alpha-beta t cell" 
                 )

list_bcells <- c("naive b cell", "memory b cell")

list_monocytes <- c("classical monocyte", "non-classical monocyte", "monocyte")

list_nkTcells <- c("mature NK T cell", "type I NK T cell")

list_DCs <- c("cd141-positive myeloid dendritic cell", "plasmacytoid dendritic cell")

list_neutrophils <- c("neutrophil", "cd24 neutrophil", "nampt neutrophil")

tabular_sapiens_blood_meta$level1_annotation <- as.character(tabular_sapiens_blood_meta$level1_annotation)

tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_tcells] <- "T-cell"
tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_monocytes] <- "monocyte"
tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_nkTcells] <- "NK T cell"
tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_DCs] <- "DC"
tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_bcells] <- "B cell"
tabular_sapiens_blood_meta$level1_annotation[tabular_sapiens_blood_meta$cell_ontology_class %in% list_neutrophils] <- "neutrophils"


tabular_sapiens_blood_meta$level1_annotation <- as.factor(tabular_sapiens_blood_meta$level1_annotation)

tabular_sapiens_blood@meta.data <- tabular_sapiens_blood_meta

# saveRDS(tabular_sapiens_blood, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/Tabular_sapiens_blood_withLevel1Annotation.rds")

# Right now the `tabular_sapiens_blood` still has Ensembl ID, need to convert to gene name to match with our data

library(org.Hs.eg.db)

tabular_sapiens_blood_count <- tabular_sapiens_blood@assays$RNA@counts

master_gene_table <- mapIds(org.Hs.eg.db, keys = row.names(tabular_sapiens_blood_count), keytype = "ENSEMBL", column="SYMBOL")
master_gene_table <- as.data.frame(master_gene_table)

tabular_sapiens_blood_editedGeneNames <- tabular_sapiens_blood

row.names(tabular_sapiens_blood_count) <- master_gene_table$master_gene_table

# Remove gene with names `NA`
tabular_sapiens_blood_count <- as.data.frame.matrix(tabular_sapiens_blood_count) %>% filter(is.na(row.names(tabular_sapiens_blood_count)) == FALSE)

which(is.na(row.names(tabular_sapiens_blood_count)) == TRUE)

# Get the list on identified gene names
identified_genes <- tabular_sapiens_blood_count[which(is.na(master_gene_table$master_gene_table) == FALSE), ]

master_gene_table_identifiedGenes <- master_gene_table %>% filter(is.na(master_gene_table) == FALSE)

tabular_sapiens_blood_editedGeneNames <- subset(tabular_sapiens_blood, features = rownames(identified_genes))

identical(row.names(master_gene_table_identifiedGenes), row.names(tabular_sapiens_blood_editedGeneNames))

row.names(tabular_sapiens_blood_editedGeneNames@assays$RNA@counts) <- master_gene_table_identifiedGenes$master_gene_table
row.names(tabular_sapiens_blood_editedGeneNames@assays$RNA@data) <- master_gene_table_identifiedGenes$master_gene_table

# saveRDS(tabular_sapiens_blood_editedGeneNames, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/Tabular_sapiens_blood_editedGeneNames_withLevel1Annotation.rds")

# Remove duplicated genes
which(duplicated(row.names(tabular_sapiens_blood_editedGeneNames)))

tabular_sapiens_blood_editedGeneNames_counts <- GetAssayData(tabular_sapiens_blood_editedGeneNames,
                                                             assay = "RNA")

tabular_sapiens_blood_editedGeneNames_counts <- tabular_sapiens_blood_editedGeneNames_counts[-which(duplicated(row.names(tabular_sapiens_blood_editedGeneNames))), ]


tabular_sapiens_blood_editedGeneNames_noduplicatedGenes <- subset(tabular_sapiens_blood_editedGeneNames,
                                                                  features = row.names(tabular_sapiens_blood_editedGeneNames_counts))

tabular_sapiens_blood_editedGeneNames_noduplicatedGenes[["RNA"]]@meta.features <- data.frame(row.names = rownames(tabular_sapiens_blood_editedGeneNames_noduplicatedGenes[["RNA"]]))

# saveRDS(tabular_sapiens_blood_editedGeneNames_noduplicatedGenes, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/Tabular_sapiens_blood_editedGeneNames_noduplicatedGenes_withLevel1Annotation.rds")


####### Perform blood-tissue annotation on (T-cells, monocyte, NK cell, B cell)
data_forQuery <- data
data_forQuery$compartment[data_forQuery$compartment == "Plasma cell"] <- "Immune"

data_forQuery.immune <- subset(data_forQuery, subset = compartment == "Immune")



######## NOTE: This step is run on Server due to high computational cost

data_list <- c(tabular_sapiens_blood_editedGeneNames_noduplicatedGenes, data_forQuery.immune)

wholeblood_data.anchor <- FindIntegrationAnchors(object.list = data_list, dims = 1:30, anchor.features = 3000)
wholeblood.integrated <- IntegrateData(anchorset = wholeblood_data.anchor, dims = 1:30)

wholeblood.integrated <- ScaleData(wholeblood.integrated, verbose = FALSE)
wholeblood.integrated <- RunPCA(wholeblood.integrated, npcs = 30, verbose = FALSE)
wholeblood.integrated <- RunUMAP(wholeblood.integrated, reduction = "pca", dims = 1:30, n.components = 10)
wholeblood.integrated <- FindNeighbors(wholeblood.integrated, reduction = "pca", dims = 1:30)
wholeblood.integrated <- FindClusters(wholeblood.integrated, resolution = 0.4)

# saveRDS(wholeblood.integrated, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/wholeblood_immune_integrated_forTBannotation.rds")

# Perform tissue and blood immune cells identification
wholeblood.integrated <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/data/wholeblood_immune_integrated_forTBannotation_viaGreatlakes2.rds")


# Define query and reference on common immune cell types
## PDAC side: NK cell, T-cell, Regulatory T-cell, M2 macrophage, B-cell, Monocyte, Plasma cell
## Reference side: nk cell, T-cell, macrophage, B cell, monocyte, plasma cell
list_common_immune_query <- c("NK cell", "T-cell", "Regulatory T-cell", "M2 macrophage",
                              "B-cell", "Monocyte", "Plasma cell")
list_common_immune_reference <- c("nk cell", "T-cell", "macrophage", "B cell", "monocyte", "plasma cell")

wholeblood.integrated_meta <- wholeblood.integrated@meta.data

wholeblood.integrated_meta$annotation_query_ref <- rep(NA, nrow(wholeblood.integrated_meta))

wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "NK cell"] <- "Query NK cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "T-cell"] <- "Query T-cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Regulatory T-cell"] <- "Query T-cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "M2 macrophage"] <- "Query macrophage"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "B-cell"] <- "Query B cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Monocyte"] <- "Query monocyte"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$annotation == "Plasma cell"] <- "Query plasma cell"

wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "nk cell"] <- "Ref NK cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "T-cell"] <- "Ref T-cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "macrophage"] <- "Ref macrophage"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "B cell"] <- "Ref B cell"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "monocyte"] <- "Ref monocyte"
wholeblood.integrated_meta$annotation_query_ref[wholeblood.integrated_meta$level1_annotation == "plasma cell"] <- "Ref plasma cell"


wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "NK cell"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "T-cell"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Regulatory T-cell"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "M2 macrophage"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "B-cell"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Monocyte"] <- "Query"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$annotation == "Plasma cell"] <- "Query"

wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "nk cell"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "T-cell"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "macrophage"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "B cell"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "monocyte"] <- "Ref"
wholeblood.integrated_meta$orig.ident[wholeblood.integrated_meta$level1_annotation == "plasma cell"] <- "Ref"


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

###Use only cells within range to avoid cluster error: -5.3 < UMAP 1 < 6.15; -2.14 < UMAP 2 < 9.65
temp1 <- temp1[temp1$UMAP_1 > -5.3 & temp1$UMAP_1 < 6.15 & temp1$UMAP_2 > -2.14 & temp1$UMAP_2 < 9.65, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_Tcells_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_Tcells_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_Tcells_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_Tcells_UMAP_d10_query_ref_dist.RData")

############# 
ct = "monocyte"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -5.4 < UMAP 1 < 3.0; -1.71 < UMAP 2 < 9.66
temp1 <- temp1[temp1$UMAP_1 > -5.4 & temp1$UMAP_1 < 3.0 & temp1$UMAP_2 > -1.71 & temp1$UMAP_2 < 9.66, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_monocyte_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_monocyte_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_monocyte_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Tissue"
km.res3$origin_tb[km.res3$cluster == 2] <- "Blood"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_monocyte_UMAP_d10_query_ref_dist.RData")

############# 
ct = "B cell"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -4.8 < UMAP 1 < 2.30; -1.05 < UMAP 2 < 10.70
temp1 <- temp1[temp1$UMAP_1 > -4.8 & temp1$UMAP_1 < 2.3 & temp1$UMAP_2 > -1.05 & temp1$UMAP_2 < 10.70, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_bcell_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_bcell_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_bcell_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_bcell_UMAP_d10_query_ref_dist.RData")

############# 
ct = "NK cell"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -5.01 < UMAP 1 < 6.10; -1.10 < UMAP 2 < 10.45
temp1 <- temp1[temp1$UMAP_1 > -5.10 & temp1$UMAP_1 < 6.10 & temp1$UMAP_2 > -1.10 & temp1$UMAP_2 < 10.45, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_NKcell_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()


###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_NKcell_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_NKcell_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_NKcell_UMAP_d10_query_ref_dist.RData")

############# 
ct = "plasma cell"

temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -2.45 < UMAP 1 < 3.95; -0.42 < UMAP 2 < 9.47
temp1 <- temp1[temp1$UMAP_1 > -2.45 & temp1$UMAP_1 < 3.95 & temp1$UMAP_2 > -0.42 & temp1$UMAP_2 < 9.47, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_plasmacell_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_plasmacell_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_plasmacell_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Tissue"
km.res3$origin_tb[km.res3$cluster == 2] <- "Blood"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_plasmacell_UMAP_d10_query_ref_dist.RData")

#############
ct = "macrophage"
temp1 <- wholeblood.integrated.common.immune_subset_meta_2[grep(ct, wholeblood.integrated.common.immune_subset_meta_2$annotation_query_ref), ]

###Use only cells within range to avoid cluster error: -5.18 < UMAP 1 < -3.50; -1.52 < UMAP 2 < -0.35
temp1 <- temp1[temp1$UMAP_1 > -5.18 & temp1$UMAP_1 < -3.50 & temp1$UMAP_2 > -1.52 & temp1$UMAP_2 < -0.35, ]
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

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_macrophage_UMAP_d10_query_ref_dist.pdf", width = 15, height = 6)
print(pheatmap(temp4, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, color = my.color1, breaks = my.break1))
dev.off()

###median distance for each query
temp4 <- as.data.frame(temp4)
temp4$dist_median = rowMedians(as.matrix(temp4[,c(-1)]))
temp5 <- temp4[, "dist_median", drop = F]
pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_macrophage_UMAP_d10_query_ref_dist_median.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE, color = my.color1, breaks = my.break1))
dev.off()

pdf(file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_macrophage_UMAP_d10_query_ref_dist_median_nobreak.pdf", width = 4, height = 6)
print(pheatmap(temp5, annotation_names_col = FALSE, annotation_names_row = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = paste0(ct), annotation_row = km.res3, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

##Assign tissue and blood origin based on distance calculated
km.res3 <- as.data.frame(km.res3, drop = F)
km.res3$origin_tb <- rep("Null", dim(km.res3)[1])
km.res3$origin_tb[km.res3$cluster == 1] <- "Blood"
km.res3$origin_tb[km.res3$cluster == 2] <- "Tissue"

save(km.res3, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_macrophage_UMAP_d10_query_ref_dist.RData")


##collect all data
temp_meta <- NULL
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_Tcells_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_monocyte_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_bcell_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_NKcell_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_plasmacell_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)
load("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/Wholeblood_macrophage_UMAP_d10_query_ref_dist.RData", verbose = T)
temp_meta <- rbind(temp_meta, km.res3)

save(temp_meta, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/tissue_blood_annotation_usingWholeBlood.RData")

### pull results and annotate tissue and blood origin
data <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/PDAC_TB_TuN_annotation_variant_analysis_4_addMonocytes_afterMonocyteTBannotated_020824.rds")
data_TB_annotation <- data

data_TB_annotation$compartment[data_TB_annotation$compartment == "Plasma cell"] <- "Immune"

data_TB_annotation_meta <- data_TB_annotation@meta.data

### tissue and blood origin assignment
data_TB_annotation_meta_1 <- data_TB_annotation_meta[!(row.names(data_TB_annotation_meta) %in% row.names(temp_meta)), ]
data_TB_annotation_meta_1$TB_origin_wholeblood <- rep("Undetermined", dim(data_TB_annotation_meta_1)[1])
data_TB_annotation_meta_1$TB_origin_wholeblood[data_TB_annotation_meta_1$compartment == "Stromal"| 
                                                 data_TB_annotation_meta_1$compartment == "Pancreas_compartment" |
                                                 data_TB_annotation_meta_1$compartment == "Beta cell"] <- "Tissue"

data_TB_annotation_meta_1$TB_origin_wholeblood[data_TB_annotation_meta_1$annotation == "MDSC"] <- "Tissue"
data_TB_annotation_meta_1$TB_origin_wholeblood[data_TB_annotation_meta_1$annotation == "Mast cell"] <- "Tissue"


data_TB_annotation_meta_2 <- data_TB_annotation_meta_1[, "TB_origin_wholeblood", drop = F]

temp_meta <- temp_meta[, "origin_tb", drop = F]
colnames(temp_meta) <- "TB_origin_wholeblood"

temp_meta <- rbind(temp_meta, data_TB_annotation_meta_2)

temp_meta <- temp_meta[row.names(data_TB_annotation_meta), , drop = F]

data_TB_annotation[["TB_origin_wholeblood"]] <- temp_meta$TB_origin_wholeblood

# Save TB annotated file
saveRDS(data_TB_annotation, file = "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/PDAC_TB_TuN_annotation_variant_analysis_4_usingWholeblood_021224.rds")


##############  Visualization

#### UMAP of all cells; cases vs. controls cells
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/UMAP_PDAC_TB_annotation_usingWholeblood.pdf", width = 12, height = 7)
DimPlot(data_TB_annotation, group.by = "annotation", label = TRUE)
DimPlot(data_TB_annotation, group.by = "annotation", label = FALSE)
DimPlot(data_TB_annotation, group.by = "annotation", split.by = "TB_origin_wholeblood", label = TRUE)
DimPlot(data_TB_annotation, group.by = "annotation", split.by = "TB_origin_wholeblood", label = FALSE)
dev.off()

#### UMAP of all cells; cases
data_TB_annotation.cases <- subset(data_TB_annotation, subset = group == "Cases")

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/UMAP_PDAC_TB_annotation_usingWholeblood_cases.pdf", width = 12, height = 7)
DimPlot(data_TB_annotation.cases, group.by = "annotation", label = TRUE)
DimPlot(data_TB_annotation.cases, group.by = "annotation", label = FALSE)
DimPlot(data_TB_annotation.cases, group.by = "annotation", split.by = "TB_origin_wholeblood", label = TRUE)
DimPlot(data_TB_annotation.cases, group.by = "annotation", split.by = "TB_origin_wholeblood", label = FALSE)
dev.off()

#### Plot UMAP plot of cell types in tumor tissues after REMOVING 
########## 1) blood immune cells, 2) non-mutated ductal and epithelial cell

data_TB_annotation.cases.removedArtifacts <- data_TB_annotation.cases

# Remove blood immune cells
data_TB_annotation.cases.removedArtifacts <- subset(data_TB_annotation.cases.removedArtifacts,
                                                  subset = TB_origin_wholeblood != "Blood")

# Remove non-tumor cells
## Keep other associated cells (TuN_Origin_new == "NA)
## Keep tumor cells (TuN_Origin_new == "Tumor cell)
data_TB_annotation.cases.removedArtifacts <- subset(data_TB_annotation.cases.removedArtifacts,
                                                  subset = TuN_Origin_new != "Non-tumor cell")


# Get only data from cases
data_TB_annotation.cases.removedArtifacts <- subset(data_TB_annotation.cases.removedArtifacts,
                                                        subset = group == "Cases")

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/UMAP_PDAC_TB_annotation_usingWholeblood_cases_afterRemoveArtifacts.pdf", width = 12, height = 10)
DimPlot(data_TB_annotation.cases.removedArtifacts, group.by = "annotation")
dev.off()

#### Plot UMAP plot of all immune cell types in blood and tissue immune cells in tumor tissues 
data_TB_annotation.cases.allImmune <- subset(data_TB_annotation.cases, 
                                                subset = compartment == "Immune")

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/UMAP_PDAC_TB_annotation_usingWholeblood_cases_allImmuneTB.pdf", width = 12, height = 10)

dimplot_all_immune_TB <- DimPlot(data_TB_annotation.cases.allImmune, group.by = "annotation", split.by = "TB_origin_wholeblood")
dimplot_all_immune_TB$data$TB_origin_wholeblood <- factor(x = dimplot_all_immune_TB$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_all_immune_TB

dev.off()


#### Plot UMAP plot of common immune cell types in blood and tissue immune cells in tumor tissues 
list_common_immune_TB <- c("T-cell", "Regulatory T-cell", "Monocyte", "B-cell", "NK cell", "Plasma cell", "M2 macrophage")

data_TB_annotation.cases.commonImmune <- subset(data_TB_annotation.cases, 
                                                subset = annotation %in% list_common_immune_TB)

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/UMAP_PDAC_TB_annotation_usingWholeblood_cases_commomImmuneTB.pdf", width = 12, height = 10)

dimplot_common_immune_TB <- DimPlot(data_TB_annotation.cases.commonImmune, group.by = "annotation", split.by = "TB_origin_wholeblood")
dimplot_common_immune_TB$data$TB_origin_wholeblood <- factor(x = dimplot_common_immune_TB$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_common_immune_TB

dev.off()


save.image("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/PDAC_TBannotation_usingWholeBlood_020924.RData")