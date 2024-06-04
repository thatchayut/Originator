#' Integrate
#' @description
#' An object with functions that help integrate query data with reference data
#'
#' @section Functions:
#' - \code{\link[=Integrate.do]{Integrate$do}}
#' - \code{\link[=Integrate.doHarmony]{Integrate$doHarmony}}
#' - \code{\link[=Integrate.annotateRefQuery]{Integrate$annotateRefQuery}}
#'
#' @export
Integrate <- new.env()

#' Integrate$do
#' @name Integrate.do
#' @description
#' Do integration. Use default Seurat pipeline (CCA), which is very slow. For
#' better options, use \code{\link[=Integrate.doHarmony]{Integrate$doHarmony}}.
#'
#' **This function is deprecated, do not use.**
#'
#' @param blood Blood reference data
#' @param query Query data
#'
#' @return Integrated data
Integrate$do <- function(blood, query) {
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
    Seurat::RunPCA(wholeblood.integrated, npcs = 30, verbose = FALSE)
  wholeblood.integrated <-
    Seurat::RunUMAP(
      wholeblood.integrated,
      reduction = "pca",
      dims = 1:30,
      n.components = 10
    )
  wholeblood.integrated <-
    Seurat::FindNeighbors(wholeblood.integrated,
                          reduction = "pca",
                          dims = 1:30)
  wholeblood.integrated <-
    Seurat::FindClusters(wholeblood.integrated, resolution = 0.4)

  return(wholeblood.integrated)

}

#' Integrate$doHarmony
#' @name Integrate.doHarmony
#' @description
#' Do integration using harmony
#'
#' @param blood Blood reference data
#' @param query Query data
#'
#' @return Integrated data
Integrate$doHarmony <- function(blood, query) {
  blood@meta.data$orig.ident <- "Ref"
  query@meta.data$orig.ident <- "Query"
  merged_data <- merge(blood, query)
  merged_data <- Seurat::NormalizeData(merged_data)
  merged_data <-
    Seurat::FindVariableFeatures(merged_data, nfeatures = 3000)
  merged_data <- Seurat::ScaleData(merged_data)
  merged_data <- Seurat::RunPCA(merged_data, npcs = 30)
  merged_data <- harmony::RunHarmony(merged_data, "orig.ident")
  merged_data <-
    Seurat::RunUMAP(
      merged_data,
      reduction = "harmony",
      dims = 1:30,
      n.components = 10
    )
  merged_data <-
    Seurat::FindNeighbors(merged_data, reduction = "harmony", dims = 1:30)
  merged_data <- Seurat::FindClusters(merged_data, resolution = 0.4)
  return(merged_data)
}

#' Integrate$annotateRefQuery
#' @name Integrate.annotateRefQuery
#' @description
#' Annotate ref and query identity for originator
#'
#' @param config Dataset config
#' @param data Integrated data
#'
#' @return Integrated data with ref and query cells annotated
Integrate$annotateRefQuery <- function(config, data) {
  meta <- data@meta.data
  # Init ref/query annotation as NA
  meta$annotation_query_ref <- rep(NA, nrow(meta))

  # Loop through all common cell type names
  for (cell_type in names(config$cell_types)) {
    # Get query cell types that are identical to the common cell type
    query_cell_types <- config$cell_types[[cell_type]]
    if (!is.list(query_cell_types)) {
      query_cell_types <- list(query_cell_types)
    }
    # Get ref cell type name that is identical to the common cell type
    ref_cell_type <- Values$wholeblood_immune_cell_type_map[[cell_type]]
    # Annotate query cells
    for (query_cell_type in query_cell_types) {
      meta$annotation_query_ref[meta[[config$ann_col]] == query_cell_type] <-
        paste0("Query ", cell_type)
    }
    # Annotate ref cells
    # TODO level1_annotation is hard coded for wholeblood data
    meta$annotation_query_ref[meta$level1_annotation == ref_cell_type] <-
      paste0("Ref ", cell_type)
  }

  data@meta.data <- meta

  # Get unique list of cell types
  # TODO Refactor this using config
  list_annotation_query_ref <-
    unlist(unique(data$annotation_query_ref))
  list_annotation_query_ref <-
    list_annotation_query_ref[!is.na(list_annotation_query_ref)]

  # Subset to get only cell types of interests
  data <-
    subset(data, subset = annotation_query_ref %in% list_annotation_query_ref)
  return(data)
}
