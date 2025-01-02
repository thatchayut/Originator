source("./src/Configs.R")
library(dplyr)


loadTissueBloodAnnotation <- function(config) {
  output_path <- paste0("./output/", config$name, "/")
  data_files <- list.files(path = output_path, pattern = "\\.rdata$", ignore.case = T)
  meta <- NULL
  for (file in data_files) {
    file_path <- paste0(output_path, file)
    print(paste("Loading", file_path))
    load(file_path)
    meta <- rbind(meta, km.res3)
  }
  return(meta)
}

runUmap <- function(data) {
  data <- Seurat::NormalizeData(data)
  data <- Seurat::FindVariableFeatures(data, nfeatures = 3000)
  data <- Seurat::ScaleData(data)
  data <- Seurat::RunPCA(data, npcs = 30)
  data <-
    Seurat::RunUMAP(
      data,
      reduction = "pca",
      dims = 1:30,
      n.components = 10
    )
  return(data)
}

runUmapIfNotExists <- function(data) {
  if (is.null(data@reductions$umap)) {
    data <- runUmap(data)
  }
  return(data)
}

dropNaCellTypeLevels <- function(config, data) {
  data@meta.data[[config$ann_col]] <- droplevels(data@meta.data[[config$ann_col]])
  return(data)
}

loadAnnotatedData <- function(config, data = NULL, annotation = NULL){
  if (is.null(data) && is.null(annotation)) {
    data <- readRDS(paste0("./data/", config$name, ".rds"))
    data <- runUmapIfNotExists(data)
    annotation <- loadTissueBloodAnnotation(config)
  }
  if (is.factor(data@meta.data[[config$ann_col]])) {
    levels(data@meta.data[[config$ann_col]]) <- unique(c(levels(data@meta.data[[config$ann_col]]), config$unified))
  }
  # data <- subset(data, subset = cell_type != "")
  for (i in 1:length(config$query)) {
    query_cell_type <- config$query[i]
    unified_cell_type <- config$unified[i]
    data[[config$ann_col]][data[[config$ann_col]] == query_cell_type] <- unified_cell_type
  }
  if (is.factor(data@meta.data[[config$ann_col]])) {
    data <- dropNaCellTypeLevels(config, data)
  }
  data$orig.ident <- "Tissue"
  data$orig.ident[row.names(annotation)] <- annotation$origin_tb
  return(data)
}

subsetCellTypesInConfig <- function(config, data) {
  data <- data[, which(data@meta.data[[config$ann_col]] %in% config$unified)]
  if (is.factor(data@meta.data[[config$ann_col]])) {
    data <- dropNaCellTypeLevels(config, data)
  }
  return(data)
}

plotTissueBloodUmap <- function(config, data) {
  data_sub <- subsetCellTypesInConfig(config, data)
  plots <- list(
    Seurat::DimPlot(data, reduction = "umap", group.by = config$ann_col, raster = F),
    Seurat::DimPlot(data, reduction = "umap", group.by = config$ann_col, label = T, raster = F),
    Seurat::DimPlot(data, reduction = "umap", group.by = config$ann_col, split.by = "orig.ident", raster = F),
    Seurat::DimPlot(data, reduction = "umap", group.by = config$ann_col, split.by = "orig.ident", label = T, raster = F),
    Seurat::DimPlot(data_sub, reduction = "umap", group.by = config$ann_col, split.by = "orig.ident", raster = F),
    Seurat::DimPlot(data_sub, reduction = "umap", group.by = config$ann_col, split.by = "orig.ident", label = T, raster = F)
  )
  return(plots)
}

calcCellTypeProportion <- function(data, ann_col) {
  celltype_proportions <- as.data.frame(prop.table(table(data[[ann_col]])))
  colnames(celltype_proportions) <- c("cellType", "proportion")
  celltype_proportions$x <- factor("proportions")
  return(celltype_proportions)
}

plotCellTypeProportion <- function(config, data) {
  plotHelper <- function(proportions) {
    colors = c(
      RColorBrewer::brewer.pal(name = "Set1", n = 9),
      RColorBrewer::brewer.pal(name = "Set3", n = 12),
      RColorBrewer::brewer.pal(name = "Set2", n = 8),
      RColorBrewer::brewer.pal(name = "Paired", n = 12),
      RColorBrewer::brewer.pal(name = "Dark2", n = 8),
      RColorBrewer::brewer.pal(name = "Accent", n = 8)
    )
    colors = c(colors, colors)
    p <- ggplot2::ggplot(proportions,
                         ggplot2::aes(x = x, y = proportion, fill = cellType)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Proportion of Cell Types", x = "", y = "Proportion") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::facet_wrap( ~ orig.ident, ncol = length(proportions))
    return(p)
  }
  
  plots = list()
  
  proportions_all <- calcCellTypeProportion(data, config$ann_col)
  subset_data_tissue <- subset(data, subset = orig.ident == "Tissue")
  proportions_tissue <- calcCellTypeProportion(subset_data_tissue, config$ann_col)
  subset_data_blood <- subset(data, subset = orig.ident == "Blood")
  proportions_blood <- calcCellTypeProportion(subset_data_blood, config$ann_col)
  
  proportions_list <- list(proportions_all, proportions_tissue, proportions_blood)
  combined_proportions <- do.call(rbind, proportions_list)
  proportions_category = c("All", "Tissue", "Blood")
  combined_proportions$orig.ident <- factor(rep(proportions_category, sapply(proportions_list, nrow)), levels = proportions_category)
  
  plots <- append(plots, list(plotHelper(combined_proportions)))
  
  subset_data_tissue <- subsetCellTypesInConfig(config, subset_data_tissue)
  proportions_tissue <- calcCellTypeProportion(subset_data_tissue, config$ann_col)
  subset_data_blood <- subsetCellTypesInConfig(config, subset_data_blood)
  proportions_blood <- calcCellTypeProportion(subset_data_blood, config$ann_col)
  
  proportions_list <- list(proportions_tissue, proportions_blood)
  combined_proportions <- do.call(rbind, proportions_list)
  proportions_category = c("Tissue", "Blood")
  combined_proportions$orig.ident <- factor(rep(proportions_category, sapply(proportions_list, nrow)), levels = proportions_category)
  
  plots <- append(plots, list(plotHelper(combined_proportions)))
  
  return(plots)
  
}

buildTissueBloodProportionHeatmap <- function(configs) {
  all_cell_types <-
    c("B-cell",
      "DC",
      "Macrophage",
      "Monocyte",
      "NK-cell",
      "T-cell",
      "Plasma-cell")
  blood_df <-
    as.data.frame(matrix(nrow = length(configs), ncol = length(all_cell_types)))
  tissue_df <-
    as.data.frame(matrix(nrow = length(configs), ncol = length(all_cell_types)))
  colnames(blood_df) <- all_cell_types
  colnames(tissue_df) <- all_cell_types
  
  for (i in 1:length(configs)) {
    config <- configs[[i]]
    print(config$name)
    data <-
      readRDS(paste0("./data/", config$name, "_tb_annotated.rds"))
    meta <- data@meta.data
    ann_col <- sym(config$ann_col)
    meta <- meta %>% filter(!!ann_col %in% config$unified)
    total_counts <- meta %>% count(!!ann_col)
    blood_counts <- meta %>%
      filter(orig.ident == "Blood") %>%
      count(!!ann_col) %>%
      rename(blood_count = n)
    tissue_counts <- meta %>%
      filter(orig.ident == "Tissue") %>%
      count(!!ann_col) %>%
      rename(tissue_count = n)
    merged_counts <- total_counts %>%
      left_join(blood_counts, by = config$ann_col) %>%
      left_join(tissue_counts, by = config$ann_col) %>%
      tidyr::replace_na(list(blood_count = 0, tissue_count = 0))
    proportion_df <- merged_counts %>%
      mutate(
        blood_percentage = (blood_count / (blood_count + tissue_count)) * 100,
        tissue_percentage = (tissue_count / (blood_count + tissue_count)) * 100
      )
    for (cell_type in all_cell_types) {
      if (cell_type %in% proportion_df[[config$ann_col]]) {
        blood_df[i, cell_type] <-
          proportion_df[proportion_df[[config$ann_col]] == cell_type, "blood_percentage"]
        tissue_df[i, cell_type] <-
          proportion_df[proportion_df[[config$ann_col]] == cell_type, "tissue_percentage"]
      }
    }
    row.names(blood_df)[[i]] <- config$name
    row.names(tissue_df)[[i]] <- config$name
  }
  return(list(blood = blood_df, tissue = tissue_df))
}

plotTissueBloodProportionHeatmap <- function(df, title = "", label = T) {
  df$dataset <- row.names(df)
  df <-
    tidyr::gather(df, key = "cell_type", value = "percentage",-dataset)
  p <-
    ggplot2::ggplot(df, ggplot2::aes(x = dataset, y = cell_type, fill = percentage)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, "YlOrRd"),
      na.value = "white",
      limits = c(0, 100)
    ) +
    ggplot2::labs(x = "Dataset", y = "Cell Type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ggtitle(title)
  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", percentage)),
      color = ifelse(df$percentage > 50, "white", "black"),
      size = 3
    )
  }
  return(p)
}