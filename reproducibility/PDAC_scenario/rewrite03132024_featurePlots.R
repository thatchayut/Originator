# This file is used to make feature plots of selected genes

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

############# T-cell
data.TB_annotation.cases.tcell <- subset(data.TB_annotation.cases, subset = annotation == "T-cell")

# Make UMAP comparing Tissue vs. Blood immune cells
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_tcell.pdf", width = 12, height = 10)

dimplot_tcell <- DimPlot(data.TB_annotation.cases.tcell, group.by = "annotation", split.by = "TB_origin_wholeblood")
dimplot_tcell$data$TB_origin_wholeblood <- factor(x = dimplot_tcell$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_tcell
dev.off()

# Make feature plot of 
## CCL4 
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_tcell_CCL4.pdf", width = 12, height = 10)

featureplot_tcell_CCL4 <- FeaturePlot(data.TB_annotation.cases.tcell, features = c("CCL4"))
featureplot_tcell_CCL4

featureplot_tcell_CCL4_2 <- FeaturePlot(data.TB_annotation.cases.tcell, features = c("CCL4"), split.by = "TB_origin_wholeblood", pt.size = 0.5)
featureplot_tcell_CCL4_2

dev.off()

## CCL5
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_tcell_CCL5.pdf", width = 12, height = 10)

featureplot_tcell_CCL5 <- FeaturePlot(data.TB_annotation.cases.tcell, features = c("CCL5"))
featureplot_tcell_CCL5

featureplot_tcell_CCL5_2 <- FeaturePlot(data.TB_annotation.cases.tcell, features = c("CCL5"), split.by = "TB_origin_wholeblood", pt.size = 0.5)
featureplot_tcell_CCL5_2

dev.off()

#############  Macrophage
data.TB_annotation.cases.macrophage <- subset(data.TB_annotation.cases, subset = annotation == "M2 macrophage")

# Make UMAP comparing Tissue vs. Blood immune cells
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_macrophage.pdf", width = 12, height = 10)

dimplot_macrophage <- DimPlot(data.TB_annotation.cases.macrophage, group.by = "annotation", split.by = "TB_origin_wholeblood", pt.size = 0.5)
dimplot_macrophage$data$TB_origin_wholeblood <- factor(x = dimplot_macrophage$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_macrophage
dev.off()

# Make feature plot of 
## INHBA
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_macrophage_INHBA.pdf", width = 12, height = 10)
featureplot_macrophage_INHBA <- FeaturePlot(data.TB_annotation.cases.macrophage, features = c("INHBA"))
featureplot_macrophage_INHBA

featureplot_macrophage_INHBA_2 <- FeaturePlot(data.TB_annotation.cases.macrophage, features = c("INHBA"), split.by = "TB_origin_wholeblood", pt.size = 1.0)
featureplot_macrophage_INHBA_2

dev.off()

#############  NK cell
data.TB_annotation.cases.nk <- subset(data.TB_annotation.cases, subset = annotation == "NK cell")

# Make UMAP comparing Tissue vs. Blood immune cells
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_NKcell.pdf", width = 12, height = 10)

dimplot_nk <- DimPlot(data.TB_annotation.cases.nk, group.by = "annotation", split.by = "TB_origin_wholeblood", pt.size = 0.5)
dimplot_nk$data$TB_origin_wholeblood <- factor(x = dimplot_nk$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_nk
dev.off()

# Make feature plot of 
## GZMK
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_NKcell_GZMK.pdf", width = 12, height = 10)
featureplot_nk_GZMK <- FeaturePlot(data.TB_annotation.cases.nk, features = c("GZMK"))
featureplot_nk_GZMK

featureplot_nk_GZMK_2 <- FeaturePlot(data.TB_annotation.cases.nk, features = c("GZMK"), split.by = "TB_origin_wholeblood", pt.size = 0.5)
featureplot_nk_GZMK_2

dev.off()


## GZMB
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_NKcell_GZMB.pdf", width = 12, height = 10)
featureplot_nk_GZMB <- FeaturePlot(data.TB_annotation.cases.nk, features = c("GZMB"))
featureplot_nk_GZMB

featureplot_nk_GZMB_2 <- FeaturePlot(data.TB_annotation.cases.nk, features = c("GZMB"), split.by = "TB_origin_wholeblood", pt.size = 0.5)
featureplot_nk_GZMB_2

dev.off()

## PRF1
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_NKcell_PRF1.pdf", width = 12, height = 10)
featureplot_nk_PRF1 <- FeaturePlot(data.TB_annotation.cases.nk, features = c("PRF1"))
featureplot_nk_PRF1

featureplot_nk_PRF1_2 <- FeaturePlot(data.TB_annotation.cases.nk, features = c("PRF1"), split.by = "TB_origin_wholeblood", pt.size = 1.0)
featureplot_nk_PRF1_2

dev.off()

#############  Monocyte
data.TB_annotation.cases.monocyte <- subset(data.TB_annotation.cases, subset = annotation == "Monocyte")

# Make UMAP comparing Tissue vs. Blood immune cells
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/UMAP_PDAC_TB_annotation_usingWholeblood_cases_monocyte.pdf", width = 12, height = 10)

dimplot_monocyte <- DimPlot(data.TB_annotation.cases.monocyte, group.by = "annotation", split.by = "TB_origin_wholeblood", pt.size = 1.0)
dimplot_monocyte$data$TB_origin_wholeblood <- factor(x = dimplot_monocyte$data$TB_origin_wholeblood, levels = c("Tissue", "Blood"))
dimplot_monocyte
dev.off()

# Make feature plot of 
## LILRA5
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/FeaturePlot_PDAC_TB_annotation_usingWholeblood_cases_monocyte_LILRA5.pdf", width = 12, height = 10)
featureplot_monocyte_LILRA5 <- FeaturePlot(data.TB_annotation.cases.monocyte, features = c("LILRA5"))
featureplot_monocyte_LILRA5

featureplot_monocyte_LILRA5_2 <- FeaturePlot(data.TB_annotation.cases.monocyte, features = c("LILRA5"), split.by = "TB_origin_wholeblood", pt.size = 2.0)
featureplot_monocyte_LILRA5_2

dev.off()


