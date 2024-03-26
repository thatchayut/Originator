# This file is used to make barplot of cell type proportion 

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

library(clustree)

library(tidyr)
library(RColorBrewer)
library(pheatmap)

library(CellChat)

library(ggVennDiagram)

library(SeuratDisk)

library(caret)

library(UpSetR)

setwd("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation")

`%!in%` = Negate(`%in%`)

# Read input data
data.TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/PDAC_TB_TuN_annotation_variant_analysis_4_usingWholeblood_03132024.rds")

# Prepare the tumor tissue data after removing artifacts
data.TB_annotation.removedArtifacts <- data.TB_annotation

## Remove blood immune cells
data.TB_annotation.removedArtifacts <- subset(data.TB_annotation.removedArtifacts,
                                                  subset = TB_origin_wholeblood != "Blood")


# Get only data from cases
data.TB_annotation.removedArtifacts.cases <- subset(data.TB_annotation.removedArtifacts,
                                                        subset = group == "Cases")


###### Get cell type proportion
tmp_data.TB_annotation.cases <- subset(data.TB_annotation, subset = group == "Cases")

tmp_data.TB_annotation.cases.bloodImmune <- subset(tmp_data.TB_annotation.cases, subset = TB_origin_wholeblood == "Blood")
tmp_data.TB_annotation.cases.bloodImmune <- subset(tmp_data.TB_annotation.cases.bloodImmune, 
                                                       subset = compartment %in% c("Immune"))

tmp_data.TB_annotation.cases.tissueImmune <- subset(tmp_data.TB_annotation.cases, subset = TB_origin_wholeblood == "Tissue")
tmp_data.TB_annotation.cases.tissueImmune <- subset(tmp_data.TB_annotation.cases.tissueImmune,
                                                        subset = compartment %in% c("Immune"))

tmp_data.TB_annotation.beforeRemovedArtifacts.cases <- subset(data.TB_annotation, subset = group == "Cases")
tmp_data.TB_annotation.beforeRemovedArtifacts.cases$artifactStatus <- rep("beforeRemove", nrow(tmp_data.TB_annotation.beforeRemovedArtifacts.cases@meta.data))


tmp_data.TB_annotation.removedArtifacts.cases <- data.TB_annotation.removedArtifacts.cases
tmp_data.TB_annotation.removedArtifacts.cases$artifactStatus <- rep("afterRemove", nrow(tmp_data.TB_annotation.removedArtifacts.cases@meta.data))

tmp_data.TB_annotation.cases.immune <- subset(tmp_data.TB_annotation.cases,
                                                  subset = compartment %in% c("Immune"))

# Prepare cell counts table
celltype_counts_table_ggplot <- as.data.frame(table(tmp_data.TB_annotation.cases.immune$annotation, 
                                                    tmp_data.TB_annotation.cases.immune$TB_origin_wholeblood))

celltype_counts_table_ggplot_beforeRemovedArtifacts <- as.data.frame(table(tmp_data.TB_annotation.beforeRemovedArtifacts.cases$annotation, 
                                                                           tmp_data.TB_annotation.beforeRemovedArtifacts.cases$artifactStatus))

celltype_counts_table_ggplotremovedArtifacts.cases <- as.data.frame(table(tmp_data.TB_annotation.removedArtifacts.cases$annotation, 
                                                                          tmp_data.TB_annotation.removedArtifacts.cases$artifactStatus))

# Combine cell type table
celltype_counts_table_ggplot <- rbind(celltype_counts_table_ggplot, celltype_counts_table_ggplot_beforeRemovedArtifacts)

celltype_counts_table_ggplot <- rbind(celltype_counts_table_ggplot, celltype_counts_table_ggplotremovedArtifacts.cases)


celltype_counts_table_ggplot <- celltype_counts_table_ggplot %>% filter(celltype_counts_table_ggplot$Var1 != "noise")

celltype_counts_table_ggplot$Var1 <- as.character(celltype_counts_table_ggplot$Var1)

list_colors <- c("#FFD700", "#C0C0C0", "#F4A460", "#E6E6FA", "#228B22", "#40E0D0", 
                 "#CD5C5C", "#000080", "#800000", "#DA70D6", "#6A5ACD", "#808000",
                 "#FF6347", "#4682B4", "#B8860B", "#8A2BE2")

list_ordered_labels_in_plot <- c('Ductal cell', 'Acinar cell', 'Beta cell', 'Peri-islet Schwann cell', 
                                 'Epithelial cell', 'Endothelial cell', 'Fibroblast', "Plasma cell",
                                 'MDSC', 'Mast cell', 'M2 macrophage', "Monocyte", "T-cell", "Regulatory T-cell",
                                 "B-cell", "NK cell")


pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/rewrite03132024_PDAC_TBwholeblood_celltypeProportion_overall_reordered.pdf", width = 15, height = 10)
ggplot(celltype_counts_table_ggplot, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = list_colors,
                    breaks = list_ordered_labels_in_plot) 
dev.off()

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/rewrite03132024_PDAC_TBwholeblood_celltypeProportion_overall_alphabetiaclOrdered.pdf", width = 15, height = 10)
ggplot(celltype_counts_table_ggplot, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = list_colors) 
dev.off()


