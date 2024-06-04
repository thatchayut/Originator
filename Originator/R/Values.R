#' Values
#' @description
#' An object with const values
#'
#' @section Vals:
#' - \code{\link[=Values.common_immune_cells]{Values$common_immune_cells}}
#' - \code{\link[=Values.wholeblood_immune_cell_type_map]{Values$wholeblood_immune_cell_type_map}}
#'
#' @export
Values <- new.env()

#' Values$common_immune_cells
#' @name Values.common_immune_cells
#' @description
#' Characters of common immune cell names
#'
Values$common_immune_cells <- c("B-cell",
                                "DC",
                                "Macrophage",
                                "Monocyte",
                                "NK-cell",
                                "T-cell",
                                "Plasma-cell")

#' Values$wholeblood_immune_cell_type_map
#' @name Values.wholeblood_immune_cell_type_map
#' @description
#' Map of cell type names from common immune cells to the corresponding
#' annotation in the wholeblood data
#'
Values$wholeblood_immune_cell_type_map <- list(
  `B-cell` = "B cell",
  `DC` = "DC",
  `Macrophage` = "macrophage",
  `Monocyte` = "monocyte",
  `NK-cell` = "nk cell",
  `T-cell` = "T-cell",
  `Plasma-cell` = "plasma cell"
)
