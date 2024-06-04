#' DataLoader
#' @description
#' An object with functions that help load data
#'
#' @section Functions:
#' - \code{\link[=DataLoader.loadTissueBloodAnnotation]{DataLoader$loadTissueBloodAnnotation}}
#' - \code{\link[=DataLoader.loadTissueBloodAnnotatedData]{DataLoader$loadTissueBloodAnnotatedData}}
#'
#' @export
DataLoader <- new.env()

#' DataLoader$loadTissueBloodAnnotation
#' @name DataLoader.loadTissueBloodAnnotation
#' @description
#' Load and gather tissue and blood annotation data from Originator RData output
#'
#' Originator outputs must be stored under `output/${config$name}/`
#'
#' @param config Dataset config
#'
#' @return Meta data as a data frame. The data frame has cell name as row name,
#' with two columns `cluster` and `origin_tb`. `cluster` is either 1 or 2, and
#' `origin_tb` is either "Blood" or "Tissue". Refer to `origin_tb` for the
#' annotation, do not use `cluster` to annotate blood or tissue cell.
DataLoader$loadTissueBloodAnnotation <- function(config) {
  # List all RData files in the output directory
  output_path <- paste0("./output/", config$name, "/")
  data_files <-
    list.files(path = output_path,
               pattern = "\\.rdata$",
               ignore.case = T)

  meta <- NULL
  for (file in data_files) {
    file_path <- paste0(output_path, file)
    # Load data from each RData, which will be stored into `km.res3`
    Logger$info(paste("load", file_path))
    load(file_path)
    meta <- rbind(meta, km.res3)
  }
  return(meta)
}

#' DataLoader$loadTissueBloodAnnotatedData
#' @name DataLoader.loadTissueBloodAnnotatedData
#' @description
#' Load tissue and blood annotation into a Seurat object
#'
#' @param config Dataset config
#' @param data Original data without tissue blood annotation. If not provided,
#' data will be loaded from `data/${config$name}.rds`.
#'
DataLoader$loadTissueBloodAnnotatedData <- function(config, data = NULL) {
  # Load data if not provided
  if (is.null(data)) {
    data <- readRDS(paste0("./data/", config$name, ".rds"))
    data <- Plot$Umap$runIfNotExists(data)
  }
  # Load annotation data frame
  annotation <- DataLoader$loadTissueBloodAnnotation(config)
  # TODO u
  # Update factor levels of annotation column, otherwise renaming annotation to
  # a non-existing factor will raise an exception
  levels(data@meta.data[[config$ann_col]]) <- unique(c(levels(data@meta.data[[config$ann_col]]), names(config$cell_types)))
  # TODO
  # Remove data that has empty string as cell types, otherwise `Seurat::DimPlot`
  # will refuse to plot
  # data <- subset(data, subset = cell_type != "")
  # TODO f
  # Rename cell types to common names
  for (cell_type in names(config$cell_types)) {
    query_cell_types <- config$cell_types[[cell_type]]
    if (!is.list(query_cell_types)) {
      query_cell_types <- list(query_cell_types)
    }
    for (query_cell_type in query_cell_types) {
      data[[config$ann_col]][data[[config$ann_col]] == query_cell_type] <- cell_type
    }
  }
  # Drop cell types that has zero occurrences
  data <- Utils$dropNaCellTypeLevels(config, data)
  # Init orig.ident as issue resident
  data$orig.ident <- "Tissue"
  # Update orig.ident according to annotation
  data$orig.ident[row.names(annotation)] <- annotation$origin_tb
  return(data)
}
