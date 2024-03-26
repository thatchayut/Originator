# This file is used to make UMAP plots 

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

data.TB_annotation.cases <- subset(data.TB_annotation, subset = group == "Cases")

# Prepare the tumor tissue data after removing artifacts
data.TB_annotation.removedArtifacts <- data.TB_annotation

## Remove blood immune cells
data.TB_annotation.removedArtifacts <- subset(data.TB_annotation.removedArtifacts,
                                              subset = TB_origin_wholeblood != "Blood")


# Get only data from cases
data.TB_annotation.removedArtifacts.cases <- subset(data.TB_annotation.removedArtifacts,
                                                    subset = group == "Cases")


ct_level_withArtifacts <- c("Acinar cell", "Beta cell", "Ductal cell", "Endothelial cell", "Epithelial cell",
                            "Fibroblast", "M2 macrophage", "Mast cell", "MDSC", "Monocyte",
                            "NK cell", "Peri-islet Schwann cell", "Plasma cell", "Regulatory T-cell",
                            "T-cell", "B-cell")
ct_level_withoutArtifacts <- c("Acinar cell", "Beta cell", "Ductal cell", "Endothelial cell", "Epithelial cell",
                               "Fibroblast", "M2 macrophage", "Mast cell", "MDSC", "Monocyte",
                               "NK cell", "Peri-islet Schwann cell", "Plasma cell", "Regulatory T-cell",
                               "T-cell","B-cell")

data.TB_annotation.cases.withArtifacts <- data.TB_annotation.cases

data.TB_annotation.cases.withArtifacts$annotation <- factor(x = data.TB_annotation.cases.withArtifacts$annotation,
                                                            levels = ct_level_withArtifacts)


data.TB_annotation.removedArtifacts.cases$annotation <- factor(x = data.TB_annotation.removedArtifacts.cases$annotation,
                                                               levels = ct_level_withoutArtifacts)


# Pull color list to make sure that the color scheme is the same in case some cell types are missing
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=16)


# UMAP before removing artifacts
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_beforeRemoveArtifacts.pdf", width = 12, height = 10)
DimPlot(data.TB_annotation.cases.withArtifacts, group.by = "annotation", label = TRUE, cols = color_list)
DimPlot(data.TB_annotation.cases.withArtifacts, group.by = "annotation", label = FALSE, cols = color_list)

dev.off()

# UMAP after removing artifacts
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_afterRemoveArtifacts.pdf", width = 12, height = 10)
DimPlot(data.TB_annotation.removedArtifacts.cases, group.by = "annotation", label = TRUE, cols = color_list)
DimPlot(data.TB_annotation.removedArtifacts.cases, group.by = "annotation", label = FALSE, cols = color_list)

dev.off()

data_TB_annotation.cases.allImmune <- subset(data.TB_annotation.cases, 
                                             subset = compartment == "Immune")

color_list_immune <- c(color_list[16], color_list[7], color_list[8], color_list[9], color_list[10], color_list[11],
                       color_list[13], color_list[14], color_list[15])

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_allImmuneTB.pdf", width = 12, height = 10)

dimplot_all_immune_TB_withLable <- DimPlot(data_TB_annotation.cases.allImmune, group.by = "annotation", split.by = "TB_origin_wholeblood", label = TRUE, cols = color_list_immune)
dimplot_all_immune_TB_withLable$data$TB_origin_wholeblood <- factor(x = dimplot_all_immune_TB_withLable$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_all_immune_TB_withLable

dimplot_all_immune_TB <- DimPlot(data_TB_annotation.cases.allImmune, group.by = "annotation", split.by = "TB_origin_wholeblood", cols = color_list_immune)
dimplot_all_immune_TB$data$TB_origin_wholeblood <- factor(x = dimplot_all_immune_TB$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_all_immune_TB

dev.off()



