curatedAtlasToSeurat <-
  function(mtx_path, cell_info_path, gene_name_path, source_col) {
    matrix <- Matrix::readMM(mtx_path)
    cell_info <- read.csv(cell_info_path)
    gene_names <- read.table(gene_name_path, header = F)
    stopifnot(dim(matrix)[1] == length(gene_names$V1) &&
                dim(matrix)[2] == dim(cell_info)[1])
    matrix@Dimnames[[1]] <- gene_names$V1
    matrix@Dimnames[[2]] <- cell_info$cell_name
    data <- Seurat::CreateSeuratObject(counts = matrix)
    data$source <- cell_info[[source_col]]
    data$cell_type <- cell_info$cell_type
    return(data)
  }

convertCuratedAtlas <- function() {
  dir <- "./data/Data_Kidney/Data_Krishna2021_Kidney/"
  return(curatedAtlasToSeurat(
    paste0(dir, "Exp_data_UMIcounts.mtx"),
    paste0(dir, "Cells.csv"),
    paste0(dir, "Genes.txt"),
    "type"
  ))
}