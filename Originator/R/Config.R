#' Config
#' @description
  #' A config prototype wrapper
  #'
#' @field name A character string representing the name of the dataset configuration.
#' @field uri A character string representing the URI where the dataset comes from.
#' @field ann_col A character string representing the name of the annotation column. When calling `data@meta.data[[ann_col]]` from a Seurat object, the column with annotation should be returned.
#' @field cell_types A named list containing mapping of common immune cell types to annotation in the dataset.
#' @export
Config <- R6::R6Class(
  "Config",
  public = list(
    name = NULL,
    uri = NULL,
    ann_col = NULL,
    cell_types = NULL,
    #' @param name A character string representing the name of the dataset configuration.
    #' @param uri A character string representing the URI where the dataset comes from.
    #' @param ann_col A character string representing the name of the annotation column. When calling `data@meta.data[[ann_col]]` from a Seurat object, the column with annotation should be returned.
    #' @param cell_types A named list containing mapping of common immune cell types to annotation in the dataset.
    initialize = function(name = NULL, uri = NULL, ann_col = "cell_type", cell_types = NULL) {
      self$name <- name
      self$uri <- uri
      self$ann_col <- ann_col
      if (!all(names(cell_types) %in% Values$common_immune_cells)) {
        stop("cell type map has invalid key")
      }
      self$cell_types <- cell_types
    },
    #' @description
      #' Create a new configuration object with updated attributes
      #'
    #' @param ... Additional parameters to update the configuration attributes.
    #' @return A cloned Config object with updated attributes.
    cloneWith = function(...) {
      cloned_obj <- self$clone()
      args <- list(...)
      for (arg_name in names(args)) {
        cloned_obj[[arg_name]] <- args[[arg_name]]
      }
      return(cloned_obj)
    },
    #' @description
      #' override print
      #'
    print = function() {
      s <- paste0("Config \"", self$name, "\" with ", length(self$cell_types), " cell types in column \"", self$ann_col, "\"")
      cat(s)
    }
  )
)
