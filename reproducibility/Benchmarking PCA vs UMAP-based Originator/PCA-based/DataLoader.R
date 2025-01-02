DataLoader <- R6::R6Class(
  classname = "DataLoader",
  public = list(
    initialize = function(config) {
      private$config <- config
    }
  ),
  private = list(
    config = NULL
  )
)

# DataLoader$set("private", "loadTissueBloodAnnotation", function(output_path) {
#   data_files <- list.files(path = output_path, pattern = "\\.rdata$", ignore.case = T)
#   meta <- NULL
#   for (file in data_files) {
#     file_path <- paste0(output_path, file)
#     print(paste("Loading", file_path))
#     load(file_path)
#     meta <- rbind(meta, km.res3)
#   }
#   return(meta)
# })
#
# DataLoader$set("public", "loadTissueBloodAnnotatedData", function(config, data = NULL, annotation = NULL) {
#   if (is.null(data) && is.null(annotation)) {
#     data <- readRDS(data_path)
#     data <- runUmapIfNotExists(data)
#     annotation <- loadTissueBloodAnnotation(output_path)
#   }
#   levels(data@meta.data[[query_celltype_col]]) <- unique(c(levels(data@meta.data[[query_celltype_col]]), unified_celltypes))
#   # data <- subset(data, subset = cell_type != "")
#   for (i in 1:length(query_celltypes)) {
#     query_cell_type <- query_celltypes[[i]]
#     unified_cell_type <- unified_celltypes[[i]]
#     data[[query_celltype_col]][data[[query_celltype_col]] == query_cell_type] <- unified_cell_type
#   }
#   data <- dropNaCellTypeLevels(config, data)
#   data$orig.ident <- "Tissue"
#   data$orig.ident[row.names(annotation)] <- annotation$origin_tb
#   return(data)
# })
