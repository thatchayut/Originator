#' Perform blood and tissue-resident immune cell identification using Originator
#'
#' @param data Data
#' @param celtype Cell type of interest
#' @param output_dir Output path to save all results
#' @param use_ident Column name indicating query or reference
#' @param use_celltype Column name indicating cell type
#' @param into_col Column to save classification result
#' @param plot Generates intermediate plots and save in output path
#' @param offset Offset of embedding borders
#'
#' @return Tissue blood origin annotated data
#'
#' @export
classifyTissueBlood <- function(
  data,
  celltype,
  output_dir,
  use_ident = "orig.ident",
  use_celltype = "unified_celltype",
  into_col = "origin_tb",
  plot = T,
  offset = -0.1
) {
  mkdir(output_dir)
  # Subset cell type of interest
  data <- subsetSeurat(data, "unified_celltype", celltype)

  # Check if both reference and query data are available
  idents <- unique(data@meta.data[[use_ident]])
  if (!("Ref" %in% idents) || !("Query" %in% idents)) {
    stop(paste0(
      "Only ",
      idents,
      " avaialable in the data. Make sure to have both Ref and Query."
    ))
  }

  if (plot) {
    message("Generate UMAP plot comparing cells from query and reference")
    pdf(file.path(
      output_dir,
      paste0(celltype, "_UMAP_annotation_query_ref.pdf")
    ),
    width = 15,
    height = 6)
    print(Seurat::DimPlot(data, group.by = use_celltype, split.by = use_ident))
    dev.off()
  }

  # Perform Blood and tissue-resident immune cell separation
  message("Extract UMAP embeddings 1-10")
  embed <- data@reductions[["umap"]]@cell.embeddings
  # Use only cells within range to avoid cluster error
  if (length(offset) == 1) {
    offset <- rep(offset, 4)
  } else if (length(offset) == 2) {
    offset <- c(offset[[1]], offset[[1]], offset[[2]], offset[[2]])
  } else if (length(offset) != 4) {
    stop("offset dimension should be 1, 2 or 4")
  }

  min_x1 <- min(embed[, 1]) + offset[[1]]
  max_x1 <- max(embed[, 1]) - offset[[2]]
  min_x2 <- min(embed[, 2]) + offset[[3]]
  max_x2 <- max(embed[, 2]) - offset[[4]]

  message(paste0("min_x1: ", min_x1))
  message(paste0("max_x1: ", max_x1))
  message(paste0("min_x2: ", min_x2))
  message(paste0("max_x2: ", max_x2))

  embed <- embed[embed[, 1] > min_x1 & embed[, 1] < max_x1 & embed[, 2] > min_x2 & embed[, 2] < max_x2, ]

  # Calculate pairwise distance matrix between query and reference embeddings
  query_embed <- embed[which(data@meta.data[[use_ident]] == "Query"), ]
  ref_embed <- embed[which(data@meta.data[[use_ident]] == "Ref"), ]
  dist_matrix <- sqrt(outer(rowSums(query_embed ^ 2), rowSums(ref_embed ^ 2), "+") - 2 * tcrossprod(query_embed, ref_embed))

  # Perform clustering
  message("Perform clustering")
  km <- kmeans(dist_matrix, 2, nstart = 25)
  clusters <- as.factor(km$cluster)
  clusters <- clusters[order(clusters)]
  clusters <- data.frame(cluster = clusters)

  # Assign blood vs. tissue-resident immune cells
  # Cluster with a larger distance to the reference is assigned as tissue-resident immune cells
  # Cluster with a smaller distance to the reference is assigned as blood immune cells
  cluster1_rowMedians <- matrixStats::rowMedians(km$centers)[1]
  cluster2_rowMedians <- matrixStats::rowMedians(km$centers)[2]

  # Assign tissue and blood origin based on distance calculated
  clusters[[into_col]] <- NA

  if (cluster1_rowMedians > cluster2_rowMedians) {
    clusters[which(clusters$cluster == 1), into_col] <- "Tissue"
    clusters[which(clusters$cluster == 2), into_col] <- "Blood"
  } else if (cluster1_rowMedians < cluster2_rowMedians) {
    clusters[which(clusters$cluster == 1), into_col] <- "Blood"
    clusters[which(clusters$cluster == 2), into_col] <- "Tissue"
  } else {
    stop(
      "Two clusters have the same distance to the reference data. Cannot separate blood and tissue-resident immune cells."
    )
  }

  # Get median distance for each query
  dist_median <- matrixStats::rowMedians(dist_matrix[rownames(clusters), ])
  breaks <- seq(0, 15)
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaks))
  if (plot) {
    message("Generate heatmaps of median distance plot between query and reference")
    pdf(file.path(
      output_dir,
      paste0(celltype, "_median_dist.pdf")
    ),
    width = 4,
    height = 6)
    print(
      pheatmap::pheatmap(
        dist_median,
        annotation_names_col = F,
        annotation_names_row = F,
        show_rownames = F,
        show_colnames = F,
        main = celltype,
        annotation_row = clusters,
        cluster_rows = F,
        cluster_cols = F,
        color = colors,
        breaks = breaks
      )
    )
    dev.off()
  }

  data@meta.data[rownames(clusters), into_col] <- clusters[[into_col]]
  saveRDS(data, file.path(output_dir, paste0(celltype, "_annotated.rds")))

  data
}

