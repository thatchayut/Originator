#' Utils
#' @description
#' An object with miscellaneous utils
#'
#' @section Functions:
#' - \code{\link[=Utils.dropNaCellTypeLevels]{Utils$dropNaCellTypeLevels}}
#' - \code{\link[=Utils.subsetCellTypesInConfig]{Utils$subsetCellTypesInConfig}}
#' - \code{\link[=Utils.calcCellTypeProportion]{Utils$calcCellTypeProportion}}
#' - \code{\link[=Utils.renameGeneNames]{Utils$renameGeneNames}}
#' - \code{\link[=Utils.clearMetaFeatures]{Utils$clearMetaFeatures}}
#'
#' @export
Utils <- new.env()

#' Utils$dropNaCellTypeLevels
#' @name Utils.dropNaCellTypeLevels
#' @description
#' Drop cell types from data that have zero observations
#'
#' @param config Dataset config
#' @param data Seurat object
#'
#' @return Updated Seurat object
Utils$dropNaCellTypeLevels <- function(config, data) {
  # TODO Modify parameters to be (data, ann_col)
  # TODO Check if meta is a factor before droplevels
  data@meta.data[[config$ann_col]] <- droplevels(data@meta.data[[config$ann_col]])
  return(data)
}

#' Utils$subsetCellTypesInConfig
#' @name Utils.subsetCellTypesInConfig
#' @description
#' Subset data to keep cell types of interests
#'
#' @param config Dataset config
#' @param data Seurat object
#'
#' @return Updated Seurat object
Utils$subsetCellTypesInConfig <- function(config, data) {
  data <- data[, which(data@meta.data[[config$ann_col]] %in% names(config$cell_types))]
  # TODO Check factor inside the function
  if (is.factor(data@meta.data[[config$ann_col]])) {
    data <- Utils$dropNaCellTypeLevels(config, data)
  }
  return(data)
}

#' Utils$calcCellTypeProportion
#' @name Utils.calcCellTypeProportion
#' @description
#' Calculate cell type proportion for simple bar plots
#'
#' @param data Seurat object
#' @param ann_col Annotation column name
#'
#' @return Percentage dataframe with 3 columns "cellType", "proportion" and "x".
#' "x" column values are all same string "proportions"
Utils$calcCellTypeProportion <- function(data, ann_col) {
  celltype_proportions <- as.data.frame(prop.table(table(data[[ann_col]])))
  colnames(celltype_proportions) <- c("cellType", "proportion")
  celltype_proportions$x <- factor("proportions")
  return(celltype_proportions)
}

#' Utils$renameGeneNames
#' @name Utils.renameGeneNames
#' @description
#' Rename gene names for a Seurat object
#'
#' Modified from \url{https://github.com/satijalab/seurat/issues/1049#issuecomment-726734005}
#'
#' @param data Seurat object
#' @param newnames List of new gene names
#'
#' @return Updated Seurat object
Utils$renameGeneNames <- function(data, newnames) {
  # FIXME RNA is hard-coded
  RNA <- data@assays$RNA

  # Check new gene names has the same length as current data
  if (nrow(RNA) == length(newnames)) {
    # Modify counts
    if (length(RNA@counts)) {
      RNA@counts@Dimnames[[1]] <- newnames
      row.names(RNA@counts) <- newnames
    }
    # Modify data
    if (length(RNA@data)) {
      RNA@data@Dimnames[[1]] <- newnames
      row.names(RNA@data) <- newnames
    }
    # Modify scale.data
    if (length(RNA@scale.data)) {
      RNA@scale.data@Dimnames[[1]] <- newnames
      row.names(RNA@scale.data) <- newnames
    }
  } else {
    stop("New gene names must be of the same length as old one")
  }
  data@assays$RNA <- RNA

  return(data)
}

#' Utils$clearMetaFeatures
#' @name Utils.clearMetaFeatures
#' @description
#' Clear meta.features for a Seurat object. Run this before integration step.
#'
#' See \url{https://github.com/satijalab/seurat/issues/2317#issuecomment-555577298}
#'
#' @param data Seurat object
#' @return Updated Seurat object
Utils$clearMetaFeatures <- function (data) {
  # FIXME RNA is hard-coded
  data@assays$RNA@meta.features <-
    data.frame(row.names = row.names(data@assays$RNA))
  return(data)
}
