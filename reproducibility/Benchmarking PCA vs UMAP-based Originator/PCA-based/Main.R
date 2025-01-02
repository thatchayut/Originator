source("/nfs/dcmb-lgarmire/yangiwen/workspace/common/Utils.R")

source("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Integrate.R")
source("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Originator.R")
source("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Configs.R")
source("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Convert.R")
source("/nfs/dcmb-lgarmire/thatchau/originator_GB_revision/src/Visualize.R")

spec <- matrix(
  c(
    "name", "n", 1, "character",
    "path", "p", 1, "character",
    "reference", "r", 1, "character",
    "output", "o", 1, "character",
    "query-celltype-col", "qcc", 1, "character",
    "ref-celltype-col", "rcc", 1, "character",
    "query-celltypes", "qc", 1, "character",
    "ref-celltypes", "rc", 1, "character",
    "unified-celltypes", "uc", 1, "character",
    "seed", "s", 1, "integer"
  ),
  byrow = T,
  ncol = 4
)

opt <- getopt::getopt(spec)

data_name <- opt[["name"]]
data_path <- opt[["path"]]
ref_path <- opt[["reference"]]
output_path <- opt[["output"]]
query_celltype_col <- opt[["query-celltype-col"]]
ref_celltype_col <- opt[["ref-celltype-col"]]
query_celltypes <- lsplit(opt[["query-celltypes"]], "\\|")
ref_celltypes <- lsplit(opt[["ref-celltypes"]], "\\|")
unified_celltypes <- lsplit(opt[["unified-celltypes"]], "\\|")
seed <- opt[["seed"]]
if (is.null(seed)) {
  seed <- 0
}
message(paste0("Set seed ", seed))
set.seed(seed)

# originator <- Originator$new(config = humanliver)

runGeneNameConvertor <- function(config) {
  print("read data")
  data <- readRDS(paste0("./data/", config$name, ".rds"))
  print("convert gene names")
  data <- convertGeneNames(data)
  data <- clearMetaFeatures(data)
  print("save converted data")
  saveRDS(data, paste0("./data/", config$name, "_converted.rds"))
  return(data)
}

runIntegrator <- function(output_path, data_path = NULL, data = NULL, ref_path = NULL, seed = 42) {
  if (is.null(data)) {
    # Read data
    message("read query data")
    data <- readRDS(data_path)
  }
  message("read ref data")
  blood <- readRDS(ref_path)

  message("integrate data")
  integrated_data <- integrateBloodDataWithHarmony(blood, data, "orig.ident", seed)

  message("save integrated data")
  saveRDS(integrated_data, file.path(output_path, "integrated.rds"))

  return(integrated_data)
}

runOriginator <- function(
  output_path,
  query_celltype_col,
  query_celltypes,
  ref_celltype_col,
  ref_celltypes,
  unified_celltypes,
  integrated_data = NULL
) {
  if (is.null(integrated_data)) {
    message("read integrated data")
    integrated_data <- readRDS(file.path(output_path, "integrated.rds"))
  }

  message("annotate ref and query")
  annotated_data <- annotateRefQuery(
    integrated_data,
    query_celltype_col,
    query_celltypes,
    ref_celltype_col,
    ref_celltypes,
    unified_celltypes
  )

  for (celltype in unique(unified_celltypes)) {
    message(celltype)
    classification <- originator(
      annotated_data,
      celltype,
      plot = T,
      output_dir = output_path
      # offset = c(-0.1, 9.5, -0.1, -0.1)
    )
    annotated_data@meta.data[rownames(classification), "origin_tb"] <- classification$origin_tb
    saveRDS(annotated_data, file.path(output_path, paste0(celltype, "_annotated.rds")))
  }
}

runLoom <- function() {
  data <-
    SeuratDisk::Connect("./data/sc-landscape-human-liver-10XV2.loom")
  print(data)
  data <- Seurat::as.Seurat(data)
  print(data)
  saveRDS(data, "./data/liver.rds")
}

runVisualizer <- function() {
  # configs <- list(liver_guilliam)
  configs <- list(originator$config)
  for (config in configs) {
    annotated_data_path <- paste0("./data/", config$name, "_tb_annotated.rds")
    print(config$name)
    data <- NULL
    if (file.exists(annotated_data_path)) {
      data <- readRDS(paste0("./data/", config$name, "_tb_annotated.rds"))
    } else {
      data <- originator$dataLoader$loadTissueBloodAnnotatedData(config)
      saveRDS(data, paste0("./data/", config$name, "_tb_annotated.rds"))
    }
    umapPlots <- plotTissueBloodUmap(config, data)
    proportionPlots <- plotCellTypeProportion(config, data)
    pdf(
      paste0(
        "./output/",
        config$name,
        "/tissue_blood_annotation_umap.pdf"
      ),
      width = 18,
      height = 9
    )
    for (plot in umapPlots) {
      print(plot)
    }
    dev.off()
    pdf(
      paste0(
        "./output/",
        config$name,
        "/tissue_blood_annotation_proportion.pdf"
      )
    )
    for (plot in proportionPlots) {
      print(plot)
    }
    dev.off()
  }
}

runTissueBloodHeatmap <- function() {
  configs <- list(hca_lung,
                  global_liver,
                  global_lung,
                  global_spleen,
                  humanliver,
                  lung_3ca_xing_normal,
                  lung_3ca_bischoff_normal,
                  kidney_3ca_krishna_normal,
                  liver_guilliam)
  heatmaps <- buildTissueBloodProportionHeatmap(configs)
  readr::write_csv(heatmaps$blood, "./data/normal_blood.csv")
  readr::write_csv(heatmaps$tissue, "./data/normal_tissue.csv")

  pdf("./output/normal_blood.pdf", width = 9, height = 9)
  print(plotTissueBloodProportionHeatmap(heatmaps$blood, "Normal Blood"))
  print(plotTissueBloodProportionHeatmap(heatmaps$blood, "Normal Blood", label = F))
  dev.off()

  pdf("./output/normal_tissue.pdf", width = 9, height = 9)
  print(plotTissueBloodProportionHeatmap(heatmaps$tissue, "Normal Tissue"))
  print(plotTissueBloodProportionHeatmap(heatmaps$tissue, "Normal Tissue", label = F))
  dev.off()

  return(0)

  configs <- list(lung_3ca_xing_luad,
                  lung_3ca_xing_nodule,
                  lung_3ca_bischoff_tumor,
                  kidney_3ca_krishna_tumor)
  heatmaps <- buildTissueBloodProportionHeatmap(configs)
  readr::write_csv(heatmaps$blood, "./data/tumor_blood.csv")
  readr::write_csv(heatmaps$tissue, "./data/tumor_tissue.csv")

  pdf("./output/tumor_blood.pdf", width = 9, height = 9)
  print(plotTissueBloodProportionHeatmap(heatmaps$blood, "Tumor Blood"))
  print(plotTissueBloodProportionHeatmap(heatmaps$blood, "Tumor Blood", label = F))
  dev.off()

  pdf("./output/tumor_tissue.pdf", width = 9, height = 9)
  print(plotTissueBloodProportionHeatmap(heatmaps$tissue, "Tumor Tissue"))
  print(plotTissueBloodProportionHeatmap(heatmaps$tissue, "Tumor Tissue", label = F))
  dev.off()
}

main <- function() {
  message(data_name)
  mkdir(output_path)
  # data <- runGeneNameConvertor(config)
  # print(data)
  data <- runIntegrator(output_path, data_path = data_path, ref_path = ref_path, seed = seed)
  runOriginator(
    output_path,
    query_celltype_col = query_celltype_col,
    query_celltypes = query_celltypes,
    ref_celltype_col = ref_celltype_col,
    ref_celltypes = ref_celltypes,
    unified_celltypes = unified_celltypes,
    integrated_data = data
  )
}


main()
# runVisualizer()
# runTissueBloodHeatmap()
# temp4 <- runOriginator(humanliver, data_i)
# temp4 <- runOriginator(liver_guilliam, data_ii)
