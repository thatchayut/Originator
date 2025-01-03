mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  } else {
    message(paste("Path", path, "already exists"))
  }
}

subsetSeurat <- function(object, key, value = NULL, operator = "equals") {
  expr <- Seurat::FetchData(object, key, clean = "none")[, 1]
  result <- NULL
  if (operator == "equals") {
    result <- object[, which(expr == value)]
  } else if (operator == "in") {
    result <- object[, which(expr %in% value)]
  } else if (operator == "notna") {
    result <- object[, which(!is.na(expr))]
  }
  return(result)
}

#' Use harmony to integrate query data with blood reference data
#'
#' @param blood Blood reference data
#' @param query Query data
#' @param into_col Column name to indicate if sample is from query or reference
#' @param seed Random seed
#'
#' @return Integrated data
#'
#' @export
integrateBloodRef <- function(blood, query, into_col = "orig.ident", seed = 42) {
  blood@meta.data[[into_col]] <- "Ref"
  query@meta.data[[into_col]] <- "Query"
  merged_data <- merge(blood, query)
  merged_data <- Seurat::NormalizeData(merged_data)
  merged_data <-
    Seurat::FindVariableFeatures(merged_data, nfeatures = 3000)
  merged_data <- Seurat::ScaleData(merged_data)
  merged_data <- Seurat::RunPCA(merged_data, npcs = 30, seed.use = seed)
  merged_data <- harmony::RunHarmony(merged_data, into_col)
  merged_data <-
    Seurat::RunUMAP(
      merged_data,
      reduction = "harmony",
      dims = 1:30,
      n.components = 10,
      seed.use = seed
    )
  merged_data
}

#' Map cell type annotation
#'
#' @param data Data
#' @param use_col Cell type annotation column
#' @param into_col Column to save mapped annotation
#'
#' @return Processed data
#'
#' @export
mapCellTypes <- function(
  data,
  use_col,
  celltype_map,
  into_col = "unified_celltype"
) {
  data@meta.data[[into_col]] <- NA
  for (from_celltype in names(celltype_map)) {
    to_celltype <- celltype_map[[from_celltype]]
    data@meta.data[which(data@meta.data[[use_col]] == from_celltype), into_col] <- to_celltype
  }
  data
}

# TODO refactor
runPopscle <- function(
  vcf_path,
  barcodes_path,
  bam_path,
  popscle_path,
  output_path = "./result"
) {
  message(paste0("PATH=", Sys.getenv("PATH")))
  cmd <- paste(
    "./exec/popscle_pnas.sh",
    shQuote(vcf_path),
    shQuote(barcodes_path),
    shQuote(bam_path),
    shQuote(popscle_path),
    shQuote(output_path),
    "> ./logs/system.log 2> ./logs/error.log"
  )
  message(paste0("bash: ", cmd))
  system(cmd)
}
