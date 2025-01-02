# Seurat@4.4.0
# BiocManager
# org.Hs.eg.db
# AnnotationDbi
# dplyr
# harmony
# hhoeflin/hdf5r
# mojaveazure/loomR
# mojaveazure/seurat-disk
# readr

# Modified from https://github.com/satijalab/seurat/issues/1049#issuecomment-726734005
renameGenesSeurat <- function(obj, newnames) {
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print(
    "Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data."
  )
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) {
      RNA@counts@Dimnames[[1]] <- newnames
      row.names(RNA@counts) <- newnames
    }
    if (length(RNA@data)) {
      RNA@data@Dimnames[[1]] <- newnames
      row.names(RNA@data) <- newnames
    }
    if (length(RNA@scale.data)) {
      RNA@scale.data@Dimnames[[1]] <- newnames
      row.names(RNA@scale.data) <- newnames
    }
  } else {
    "Unequal gene sets: nrow(RNA) != nrow(newnames)"
  }
  obj@assays$RNA <- RNA

  return(obj)
}

convertGeneNames <- function (data) {
  counts <- data@assays$RNA@counts

  # Convert ensembl to symbol
  print("convert ensembl to symbol")
  master_gene_table <-
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = row.names(counts),
      keytype = "ENSEMBL",
      column = "SYMBOL"
    )
  master_gene_table <- as.data.frame(master_gene_table)

  # Update gene names
  print("update gene names")
  row.names(counts) <- master_gene_table$master_gene_table
  data <-
    renameGenesSeurat(data, master_gene_table$master_gene_table)

  # Remove gene with names `NA`
  print("remove genes with names NA")
  counts_row_names_is_na <- is.na(row.names(counts))
  master_gene_table <-
    master_gene_table[which(!counts_row_names_is_na), , drop = FALSE]
  counts <- counts[!counts_row_names_is_na, ]

  # assert which(is.na(row.names(counts)) == TRUE)

  # Get the list on identified gene names
  print("get the list on identified gene names")
  identified_genes <-
    counts[which(!is.na(master_gene_table$master_gene_table)), ]

  # Remove duplicates
  print("remove duplicates")
  # Remove duplicates from data
  data <- data[!duplicated(row.names(data)), ]

  # Remove duplicates from identified_genes
  identified_genes <-
    identified_genes[!duplicated(row.names(identified_genes)), ]

  # Check the dimensions of the unique objects
  # assert

  print("creating renamed object")
  data <- subset(data, features = row.names(identified_genes))

  return(data)

}

extractImmune <- function (data) {
  return(subset(data, subset = ann_level_1 == "Immune" & ann_level_4 != "None"))
}

#' @link https://github.com/satijalab/seurat/issues/2317#issuecomment-555577298
clearMetaFeatures <- function (data) {
  data@assays$RNA@meta.features <-
    data.frame(row.names = row.names(data@assays$RNA))
  return(data)
}

integrateBloodData <- function(blood, query) {
  data_list <- c(blood, query)

  wholeblood.anchor <-
    Seurat::FindIntegrationAnchors(
      object.list = data_list,
      dims = 1:30,
      anchor.features = 3000
    )
  wholeblood.integrated <-
    Seurat::IntegrateData(anchorset = wholeblood.anchor, dims = 1:30)

  wholeblood.integrated <-
    Seurat::ScaleData(wholeblood.integrated, verbose = FALSE)
  wholeblood.integrated <-
    Seurat::RunPCA(wholeblood.integrated, npcs = 10, verbose = FALSE)
  wholeblood.integrated <-
    Seurat::RunUMAP(
      wholeblood.integrated,
      reduction = "pca",
      dims = 1:10,
      n.components = 10
    )
  wholeblood.integrated <-
    Seurat::FindNeighbors(wholeblood.integrated,
                          reduction = "pca",
                          dims = 1:10)
  wholeblood.integrated <-
    Seurat::FindClusters(wholeblood.integrated, resolution = 0.4)

  return(wholeblood.integrated)

}

integrateBloodDataWithHarmony <- function(blood, query, batch_col = "batch", seed = 42) {
  blood@meta.data$orig.ident <- "Ref"
  query@meta.data$orig.ident <- "Query"
  merged_data <- merge(blood, query)
  merged_data <- Seurat::NormalizeData(merged_data)
  merged_data <-
    Seurat::FindVariableFeatures(merged_data, nfeatures = 3000)
  merged_data <- Seurat::ScaleData(merged_data)
  merged_data <- Seurat::RunPCA(merged_data, npcs = 10, seed.use = seed)
  merged_data <- harmony::RunHarmony(merged_data, batch_col)
  merged_data <-
    Seurat::RunUMAP(
      merged_data,
      reduction = "harmony",
      dims = 1:10,
      n.components = 10,
      seed.use = seed
    )
  # merged_data <-
  #   Seurat::FindNeighbors(merged_data, reduction = "harmony", dims = 1:30)
  # merged_data <- Seurat::FindClusters(merged_data, resolution = 0.4)
  return(merged_data)
}


annotateRefQuery <- function(
  data,
  query_celltype_col,
  query_celltypes,
  ref_celltype_col,
  ref_celltypes,
  unified_celltypes
) {
  meta <- data@meta.data

  for (i in 1:length(unified_celltypes)) {
    query_celltype <- query_celltypes[[i]]
    ref_celltype <- ref_celltypes[[i]]
    unified_celltype <- unified_celltypes[[i]]

    meta[which(meta[[query_celltype_col]] == query_celltype), "annotation_query_ref"] <- paste0("Query ", unified_celltype)
    meta[which(meta[[ref_celltype_col]] == ref_celltype), "annotation_query_ref"] <- paste0("Ref ", unified_celltype)
    meta[which(meta[[query_celltype_col]] == query_celltype), "unified_celltype"] <- unified_celltype
    meta[which(meta[[ref_celltype_col]] == ref_celltype), "unified_celltype"] <- unified_celltype
  }

  data@meta.data <- meta
  # Subset to get only cell types of interests
  data <- subsetSeurat(data, "unified_celltype", unified_celltypes, operator = "in")
  return(data)
}
