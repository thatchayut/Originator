# This file is used to plot a box plot of F1 score comparing Originator PCA-based and
#   UMAP-based in TB annotation of cleaned ccRCC & paired PBMC
# F1 score is processed by: "/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/observe_TB_annotation_krishna_et_al_2021_seurat_vs_originator.R"

library(Seurat)
library(dplyr)
library(ggsignif)
library(ggplot2)

setwd("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace")

base_output_path <- "E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/"

# Read data
df_F1_originator_pca <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_wholeblood_ref_originator_PCAbased.rds")
df_F1_originator_pca <- subset(df_F1_originator_pca, subset = cell_type != "mono")


df_F1_originator_umap <- readRDS("E:/Tays_workspace/Originator_manuscript/GB_revision_workspace/results/F1_TBannotation_cleanedccRCC_wholeblood_ref_originator_UMAPbased.rds")
# Remove monocyte as we will not show in the manuscript
df_F1_originator_umap <- subset(df_F1_originator_umap, subset = cell_type != "mono")

# combind F1 dataframes for visualization
## Add a column for either seurat-based approach or originator-based approach
df_F1_originator_pca$Embeddings <- "PCA"
df_F1_originator_umap$Embeddings <- "UMAP"

## Combine dataframe
df_F1_combined <- rbind(df_F1_originator_umap, df_F1_originator_pca)

df_F1_combined$Embeddings <- factor(df_F1_combined$Embeddings, levels=c("PCA", "UMAP"))

for (celltype in unique(df_F1_combined$cell_type)){
  # Visualize for each cell type one at a time
  pdf(paste0(base_output_path, "boxplot_F1_PCA_vs_UMAP_originator_withSignificant_", celltype, ".pdf"), width = 15, height = 12)
  
  tmp_df_F1 <- subset(df_F1_combined, subset = cell_type == celltype)
  
  ggplot2::ggplot(tmp_df_F1, ggplot2::aes(x = cell_type, y = score, fill = Embeddings)) +
    ggplot2::geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("skyblue", "lightgreen", "gray40")) +
    geom_point(aes(color= Embeddings, shape=Embeddings), alpha = 0.9, 
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
      aes(group=Embeddings),
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
}

###### Reformat to get different approaches in different x labels
ggplot2::ggplot(df_F1_combined, ggplot2::aes(x = cell_type, y = score, fill = Embeddings)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  ggplot2::stat_summary(
    fun = mean,
    geom = "text",
    ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
    vjust = -3,
    color = "black",
    size = 5
  )

##### CD4
df_F1_originator_pca_cd4 <- subset(df_F1_originator_pca, subset = cell_type == "cd4")
df_F1_originator_umap_cd4 <- subset(df_F1_originator_umap, subset = cell_type == "cd4")
df_F1_originator_cd4 <- rbind(df_F1_originator_pca_cd4, df_F1_originator_umap_cd4)

##### CD8
df_F1_originator_pca_cd8 <- subset(df_F1_originator_pca, subset = cell_type == "cd8")
df_F1_originator_umap_cd8 <- subset(df_F1_originator_umap, subset = cell_type == "cd8")
df_F1_originator_cd8 <- rbind(df_F1_originator_pca_cd8, df_F1_originator_umap_cd8)

##### NK
df_F1_originator_pca_nk <- subset(df_F1_originator_pca, subset = cell_type == "nk")
df_F1_originator_umap_nk <- subset(df_F1_originator_umap, subset = cell_type == "nk")
df_F1_originator_nk <- rbind(df_F1_originator_pca_nk, df_F1_originator_umap_nk)

##### B-cell
df_F1_originator_pca_b <- subset(df_F1_originator_pca, subset = cell_type == "b")
df_F1_originator_umap_b <- subset(df_F1_originator_umap, subset = cell_type == "b")
df_F1_originator_b <- rbind(df_F1_originator_pca_b, df_F1_originator_umap_b)

pdf(paste0(base_output_path, "F1_boxplot_ccRCC_krishna_PCA_vs_UMAP_CD4_4digits.pdf"), width = 8, height = 15)
#### With overlaying text
ggplot2::ggplot(df_F1_originator_cd4, ggplot2::aes(x = Embeddings, y = score)) +
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
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )

#### Without overlaying text
ggplot2::ggplot(df_F1_originator_cd4, ggplot2::aes(x = Embeddings, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )
dev.off()

pdf(paste0(base_output_path, "F1_boxplot_ccRCC_krishna_PCA_vs_UMAP_CD8_4digits.pdf"), width = 8, height = 15)
#### With overlaying text
ggplot2::ggplot(df_F1_originator_cd8, ggplot2::aes(x = Embeddings, y = score)) +
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
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )

#### Without overlaying text
ggplot2::ggplot(df_F1_originator_cd8, ggplot2::aes(x = Embeddings, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )
dev.off()

pdf(paste0(base_output_path, "F1_boxplot_ccRCC_krishna_PCA_vs_UMAP_nk_4digits.pdf"), width = 8, height = 15)
#### With overlaying text
ggplot2::ggplot(df_F1_originator_nk, ggplot2::aes(x = Embeddings, y = score)) +
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
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )

#### Without overlaying text
ggplot2::ggplot(df_F1_originator_nk, ggplot2::aes(x = Embeddings, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )
dev.off()


pdf(paste0(base_output_path, "F1_boxplot_ccRCC_krishna_PCA_vs_UMAP_b_4digits.pdf"), width = 8, height = 15)
#### With overlaying text
ggplot2::ggplot(df_F1_originator_b, ggplot2::aes(x = Embeddings, y = score)) +
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
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )

#### Without overlaying text
ggplot2::ggplot(df_F1_originator_b, ggplot2::aes(x = Embeddings, y = score)) +
  ggplot2::geom_boxplot(color = "skyblue") +
  ggplot2::geom_jitter(width = 0.2, color = "red") +
  stat_summary(
    aes(group=Embeddings),
    position=position_dodge(0.75),
    fun=mean,
    geom='point', color='darkred', shape=18, size=5,
    show.legend = FALSE
  ) +
  ggplot2::ylim(0.7, 1) +
  ggplot2::labs(y = "F1 score") +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 24)) +
  geom_signif(comparisons = list(c("PCA", "UMAP")), 
              map_signif_level=TRUE, textsize=10) +
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=30,face="bold"),
        axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.title.y.left = element_text(size = 40),
        legend.text=element_text(size=25),
        legend.title=element_text(size=10),
        axis.line = element_line(size=2)
  )
dev.off()

