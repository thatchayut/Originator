source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

# Test Originator pipeline
data_name <- "kidney_krishna_nk"
data_path <- "/nfs/dcmb-lgarmire/yangiwen/workspace/originator.bak/data/krishna/3ca/complete_response_all.rds"
ref_path <- "/nfs/dcmb-lgarmire/shared/public/artificially_mixing_blood_tissue/wholeblood_reference/Tabular_sapiens_blood_editedGeneNames_noduplicatedGenes_withLevel1Annotation.rds"
output_path <- "/nfs/dcmb-lgarmire/yangiwen/workspace/originator.bak/output/kidney_krishna_nk_pca"
query_celltype_col <- "cell_type"
ref_celltype_col <- "level1_annotation"
query_celltypes <- lsplit("NK_cell", "\\|")
ref_celltypes <- lsplit("nk cell", "\\|")
unified_celltypes <- lsplit("NK-cell", "\\|")

# Visualization
seeds <- c(0, 42, 64, 123, 894)
celltypes <- c("cd4", "cd8", "nk")
unified_celltypes <- c("CD4", "CD8", "NK-cell")
calcF1Score <- function(truth, pred) {
  mat <- table(truth, pred)
  tp <- mat[1, 1]
  fp <- mat[2, 1]
  fn <- mat[1, 2]
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  2 * (precision * recall) / (precision + recall)
}
f1_df <- data.frame()
for (i in 1:length(celltypes)) {
  f1_scores <- list()
  celltype <- celltypes[[i]]
  unified_celltype <- unified_celltypes[[i]]
  message(celltype)
  for (seed in seeds) {
    data <- readRDS(file.path("./output/kidney_krishna_seeds_pca_10", celltype, seed, paste0(unified_celltype, "_annotated.rds")))
    data <- subsetSeurat(data, "unified_celltype", unified_celltype)
    data@meta.data[which(data$type == "Normal"), "tissue"] <- "kidney"
    data@meta.data[which(data$type == "Tumor"), "tissue"] <- "kidney"
    data@meta.data[which(data$type == "PBMC"), "tissue"] <- "blood"
    data@meta.data[which(data$orig.ident == "Ref"), "tissue"] <- "blood"
    f1_score <- calcF1Score(data$tissue, data$origin_tb)
    f1_scores <- append(f1_scores, f1_score)
    # data_query <- subsetSeurat(data, "orig.ident", "Query")
    # data_ref <- subsetSeurat(data, "orig.ident", "Ref")
    # png(file.path("./output/kidney_krishna_seeds_pca", celltype, seed, paste0(unified_celltype, "_annotated.png")), width = 1024, height = 512)
    # # data_query <- data_query[, data_query@reductions$umap@cell.embeddings[, 1] > 1]
    # print(
    #   Seurat::DimPlot(data_query, group.by = "tissue") +
    #   Seurat::DimPlot(data_query, group.by = "origin_tb")
    # )
    # dev.off()
  }
  f1_scores <- unlist(f1_scores)
  df <- data.frame(cell_type = celltype, score = f1_scores)
  f1_df <- rbind(f1_df, df)
}
f1_df$iteration <- rep(1:5, 3)
f1_df$reduction <- "pca"
pdf("./output/kidney_krishna_seeds/f1_boxplot.pdf", width = 8, height = 6)
print(
  ggplot2::ggplot(f1_df, ggplot2::aes(x = cell_type, y = score)) +
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
)
dev.off()
pdf("./output/kidney_krishna_seeds/f1_heatmap.pdf", width = 6, height = 8)
print(
  ggplot2::ggplot(f1_df, ggplot2::aes(x = cell_type, y = iteration, fill = score)) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::coord_equal() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", score)), color = "white", size = 5) +
    # ggplot2::scale_fill_continuous(limits = c(0.5, 1)) +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", limits = c(0.5, 1), direction = 1) +
    ggplot2::ggtitle("F1 score") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 24))
)
dev.off()

pdf("./output/kidney_krishna_seeds/f1_boxplot_pca_nk.pdf", width = 8, height = 6)
print(
  ggplot2::ggplot(f1_df[f1_df$cell_type == "nk", ], ggplot2::aes(x = reduction, y = score)) +
    ggplot2::geom_boxplot(color = "skyblue") +
    ggplot2::geom_jitter(width = 0.2, color = "red") +
    ggplot2::stat_summary(
      fun = mean,
      geom = "text",
      ggplot2::aes(label = sprintf("Mean: %.2f", ggplot2::after_stat(y))),
      vjust = 2,
      color = "black",
      size = 5
    ) +
    ggplot2::stat_summary(
      fun.data = function(x)
        data.frame(y = mean(x), label = sprintf("SD: %.5f", sd(x))),
      geom = "text",
      vjust = 3.5,
      color = "black",
      size = 5
    ) +
    ggsignif::geom_signif(
      comparisons = list(c("umap", "pca")),
      map_signif_level = T,
      margin_top = 1,
      tip_length = c(1.5, 0.8),
      textsize = 13
    ) +
    ggplot2::ylim(0.5, 1) +
    ggplot2::labs(y = "F1 score") +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = 24))
)
dev.off()

# Test
blood_cells <- subsetSeurat(data, "tissue", "blood")
blood_cells <- Seurat::NormalizeData(blood_cells)
blood_cells <-
  Seurat::FindVariableFeatures(blood_cells, nfeatures = 3000)
blood_cells <- Seurat::ScaleData(blood_cells)
blood_cells <- Seurat::RunPCA(blood_cells, npcs = 30)
blood_cells <- harmony::RunHarmony(blood_cells, "orig.ident")
blood_cells <-
  Seurat::RunUMAP(
    blood_cells,
    reduction = "harmony",
    dims = 1:30,
    n.components = 10
  )
Seurat::DimPlot(blood_cells, group.by = "orig.ident")

blood_query <- subsetSeurat(data, "tissue", "lung")
blood_query <- subsetSeurat(blood_query, "orig.ident", "Query")
blood$orig.ident <- "Ref"
blood_allcelltypes <- merge(blood_query, blood)
blood_allcelltypes <- Seurat::NormalizeData(blood_allcelltypes)
blood_allcelltypes <-
  Seurat::FindVariableFeatures(blood_allcelltypes, nfeatures = 3000)
blood_allcelltypes <- Seurat::ScaleData(blood_allcelltypes)
blood_allcelltypes <- Seurat::RunPCA(blood_allcelltypes, npcs = 30)
blood_allcelltypes <- harmony::RunHarmony(blood_allcelltypes, "orig.ident")
blood_allcelltypes <-
  Seurat::RunUMAP(
    blood_allcelltypes,
    reduction = "harmony",
    dims = 1:30,
    n.components = 10
  )
Seurat::DimPlot(blood_allcelltypes, group.by = "level1_annotation", split.by = "orig.ident", label = T)

