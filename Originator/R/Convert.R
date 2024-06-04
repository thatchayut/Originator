#' Convert
#' @description
#' An object with functions that help convert between data
#'
#' @section Functions:
  #' - \code{\link[=Convert.ensemblToSymbol]{Convert$ensemblToSymbol}}
  #' - \code{\link[=Convert.curatedAtlasToSeurat]{Convert$curatedAtlasToSeurat}}
#'
#' @export
Convert <- new.env()

#' Convert$ensemblToSymbol
#' @name Convert.ensemblToSymbol
#' @description
#' Convert gene names from ensembl id to symbols
#'
#' @param data Seurat data
#'
#' @return Converted data
Convert$ensemblToSymbol <- function (data) {
  counts <- data@assays$RNA@counts

  # Convert ensembl to symbol
  Logger$info("convert ensembl to symbol")
  master_gene_table <-
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = row.names(counts),
      keytype = "ENSEMBL",
      column = "SYMBOL"
    )
  master_gene_table <- as.data.frame(master_gene_table)

  # Update gene names
  Logger$info("update gene names")
  row.names(counts) <- master_gene_table$master_gene_table
  data <- Utils$renameGeneNames(data, master_gene_table$master_gene_table)

  # Remove gene with names `NA`
  Logger$info("remove genes with names NA")
  counts_row_names_is_na <- is.na(row.names(counts))
  master_gene_table <-
    master_gene_table[which(!counts_row_names_is_na), , drop = FALSE]
  counts <- counts[!counts_row_names_is_na, ]

  # Get the list on identified gene names
  Logger$info("get the list on identified gene names")
  identified_genes <-
    counts[which(!is.na(master_gene_table$master_gene_table)), ]

  # Remove duplicates
  Logger$info("remove duplicates")
  # Remove duplicates from data
  data <- data[!duplicated(row.names(data)), ]

  # Remove duplicates from identified_genes
  identified_genes <-
    identified_genes[!duplicated(row.names(identified_genes)), ]

  # Create renamed seurat object
  Logger$info("create renamed seurat object")
  data <- subset(data, features = row.names(identified_genes))

  return(data)

}

#' Convert$curatedAtlasToSeurat
#' @name Convert.curatedAtlasToSeurat
#' @description
#' Read and convert data downloaded from
#' [3ca](https://www.weizmann.ac.il/sites/3CA/) to a Seurat object
#'
#' @param mtx_path path to the scRNA-seq expression matrix, typically named as
#' `*.mtx`
#' @param cell_info_path path to the cell information and meta data csv file,
#' typically named as `Cells.csv`
#' @param gene_name_path path to the gene name file, typically named as
#' `Genes.txt`
#' @param source_col column name in cell info csv indicating the source of the
#' cell, either from tumor or normal tissues
#'
#' @usage
#' dir <- "./data/kidney/"
#' Convert$curatedAtlasToSeurat(
#'   paste0(dir, "Exp_data_UMIcounts.mtx"),
#'   paste0(dir, "Cells.csv"),
#'   paste0(dir, "Genes.txt"),
#'   "type"
#' )
#'
#' @return Seurat data
Convert$curatedAtlasToSeurat <-
  function(mtx_path, cell_info_path, gene_name_path, source_col) {
    # Read as matrix
    matrix <- Matrix::readMM(mtx_path)
    # Read the cell information from a CSV file
    cell_info <- read.csv(cell_info_path)
    # Read gene names from a text file as a data frame
    gene_names <- read.table(gene_name_path, header = F)
    # Check dimensions
    if (dim(matrix)[1] != length(gene_names$V1) ||
        dim(matrix)[2] != dim(cell_info)[1]) {
      stop("Dimensions of the matrix do not match the provided gene and cell information.")
    }
    # Assign gene names to the rows and cell names to the columns of the matrix
    matrix@Dimnames[[1]] <- gene_names$V1
    matrix@Dimnames[[2]] <- cell_info$cell_name
    # Create Seurat object
    data <- Seurat::CreateSeuratObject(counts = matrix)
    data$source <- cell_info[[source_col]]
    data$cell_type <- cell_info$cell_type

    return(data)
  }
