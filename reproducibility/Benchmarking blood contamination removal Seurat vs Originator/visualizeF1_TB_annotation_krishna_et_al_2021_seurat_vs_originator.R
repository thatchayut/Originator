# This file is used to plot a box plot of F1 score comparing Originator-based and
#   Seurat-based in  TB annotation of cleaned ccRCC & paired PBMC
# F1 score is processed by: "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/observe_TB_annotation_krishna_et_al_2021_seurat_vs_originator.R"

library(Seurat)
library(dplyr)
library(ggsignif)
library(ggplot2)

setwd("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace")

base_output_path <- "E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/"

# Read data
df_F1_originator <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_wholeblood_ref.rds")

df_F1_seurat <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_seuratInferred_PBMC_ref.rds")

# combind F1 dataframes for visualization
## Add NA to monocytes of seurat-based approach as they are removed during blood contamination removal
df_dummy_mono <- data.frame(cell_type = rep("mono", 5),
                            score = rep(0, 5))
df_F1_seurat_with_mono <- rbind(df_F1_seurat, df_dummy_mono)

## Add a column for either seurat-based approach or originator-based approach
df_F1_originator$approach <- "Originator-based"
df_F1_seurat_with_mono$approach <- "Seurat-based"

## Combine dataframe
df_F1_combined <- rbind(df_F1_originator, df_F1_seurat_with_mono)

df_F1_combined$cell_type <- factor(df_F1_combined$cell_type, levels=unique(df_F1_combined$cell_type))
df_F1_combined$approach <- factor(df_F1_combined$approach, levels=c("Originator-based", "Seurat-based"))

####### With text
pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originator_separate.pdf"), width = 15, height = 12)

####### With text
### Originator only
ggplot2::ggplot(df_F1_originator, ggplot2::aes(x = cell_type, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = -3,
    color = "black",
    size = 5
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = -1.5,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24))


### Seurat-only
ggplot2::ggplot(df_F1_seurat, ggplot2::aes(x = cell_type, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = -3,
    color = "black",
    size = 5
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = -1.5,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24))

### Use this
ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = 1.5,
    color = "black",
    size = 5
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = 2.8,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24)) 

dev.off()

pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originator_combine.pdf"), width = 10, height = 8)
### Use this
ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("Originator-based", "Seurat-based")), 
              map_signif_level=TRUE)

dev.off()

##### NO monocyte
df_F1_combined_noMono <- subset(df_F1_combined, subset = cell_type != "mono")
df_F1_combined_noMono$cell_type <- factor(df_F1_combined_noMono$cell_type, levels=c("b", "cd4", "cd8", "nk"))


pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originator_combine_noMonocyte.pdf"), width = 15, height = 12)
### Use this
ggplot2::ggplot(df_F1_combined_noMono, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 40),
                 panel.grid.minor = element_line(size = 0.5), panel.grid.major = element_line(size = 1))

### With text
ggplot2::ggplot(df_F1_combined_noMono, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = 1.5,
    color = "black",
    size = 5
  ) +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = 2.8,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24))


dev.off()

pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originator_combine_noMono_withSignificant.pdf"), width = 15, height = 12)

ggplot2::ggplot(df_F1_combined_noMono, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = 1.5,
    color = "black",
    size = 5
  ) +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
    ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = 2.8,
    color = "black",
    size = 5
  ) +
  ggpubr::stat_compare_means(
    method = "t.test", # Or use "wilcox.test" or other methods
    aes(label = ..p.signif..), # Use stars for significance
    label.y = 0.95, # Adjust Y-position for the labels
    label.x = 1.5 # Adjust X-position (optional)
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24)) 

dev.off()

##### NO monocyte: New Seurat UMAP

# Read data
df_F1_originator_umap_new <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_wholeblood_ref_originator_UMAPbased.rds")

df_F1_seurat <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_seuratInferred_PBMC_ref.rds")

## Add a column for either seurat-based approach or originator-based approach
df_F1_originator_umap_new$approach <- "Originator-based"
df_F1_seurat$approach <- "Seurat-based"

## Combine dataframe
df_F1_combined <- rbind(df_F1_originator_umap_new, df_F1_seurat)

list_preferred_celltype <- c("b", "cd4", "cd8", "nk")

df_F1_combined$cell_type <- factor(df_F1_combined$cell_type, levels=list_preferred_celltype)
df_F1_combined$approach <- factor(df_F1_combined$approach, levels=c("Originator-based", "Seurat-based"))

pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originatorNewUMAP_combine.pdf"), width = 15, height = 12)
### Use this
ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 40),
                 panel.grid.minor = element_line(size = 0.5), panel.grid.major = element_line(size = 1))

### With text
ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = 1.5,
    color = "black",
    size = 5
  ) +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = 2.8,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24))


dev.off()


pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originatorNewUMAP_separate.pdf"), width = 15, height = 12)

####### With text
### Originator only
ggplot2::ggplot(df_F1_originator_umap_new, ggplot2::aes(x = cell_type, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.4f", ggplot2::after_stat(y))),
    vjust = -3,
    color = "black",
    size = 5
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = -1.5,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24))


### Seurat-only
ggplot2::ggplot(df_F1_seurat, ggplot2::aes(x = cell_type, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.4f", ggplot2::after_stat(y))),
    vjust = -3,
    color = "black",
    size = 5
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = -1.5,
    color = "black",
    size = 5
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24))

dev.off()

pdf(paste0(base_output_path, "boxplot_F1_seurat_vs_originatorNewUMAP_combine_withSignificant.pdf"), width = 15, height = 12)

ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = approach)) +
  ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
  geom_point(aes(color= approach, shape=approach), alpha = 0.9, 
             position = position_jitterdodge(jitter.width = 0.5))  +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = 1.5,
    color = "black",
    size = 5
  ) +
  stat_summary(
    aes(group=approach),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::stat_summary(
    fun.data = function(x)
      data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
    geom = "text",
    vjust = 2.8,
    color = "black",
    size = 5
  ) +
  ggpubr::stat_compare_means(
    method = "t.test", # Or use "wilcox.test" or other methods
    aes(label = ..p.signif..), # Use stars for significance
    label.y = 0.95, # Adjust Y-position for the labels
    label.x = 1.5 # Adjust X-position (optional)
  ) +
  ggplot2::ylim(0.5, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(text = ggplot2::element_text(size = 24)) 

dev.off()

