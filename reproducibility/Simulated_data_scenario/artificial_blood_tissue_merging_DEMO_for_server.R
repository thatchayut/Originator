# Note: This code is implemented for using on "SLURM" server

library(dplyr)
library(Seurat)
library(patchwork)

# Read processed files for merging
data.for.merging.M1 <- readRDS("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/data/for_merging/data_for_merging_M1.rds")
data.for.merging.M2 <- readRDS("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/data/for_merging/data_for_merging_M2.rds")
data.for.merging.pbmc8k <- readRDS("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/data/for_merging/data_for_merging_pbmc8k.rds")
data.for.merging.TI_NK.1 <- readRDS("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/data/for_merging/data_for_merging_TI_NK_1.rds")
data.for.merging.TI_NK.2 <- readRDS("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/data/for_merging/data_for_merging_TI_NK_2.rds")

data.merged.tissue.blood.residents.demo <- merge(data.for.merging.M1, y = c(data.for.merging.pbmc8k,
                                                                            data.for.merging.TI_NK.1
                                                ),
                                                add.cell.ids = c("M1", "pbmc8k",
                                                                 "TI_NK_1"), 
                                                project = "merged_DEMO")

data.merged.tissue.blood.residents.demo <- NormalizeData(data.merged.tissue.blood.residents.demo, normalization.method = "LogNormalize", scale.factor = 10000)

data.merged.tissue.blood.residents.demo <- FindVariableFeatures(data.merged.tissue.blood.residents.demo, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(data.merged.tissue.blood.residents.demo)
data.merged.tissue.blood.residents.demo <- ScaleData(data.merged.tissue.blood.residents.demo, features = all.genes)

data.merged.tissue.blood.residents.demo <- RunPCA(data.merged.tissue.blood.residents.demo, npcs = 30, verbose = FALSE)
data.merged.tissue.blood.residents.demo <- RunUMAP(data.merged.tissue.blood.residents.demo, reduction = "pca", dims = 1:30)
data.merged.tissue.blood.residents.demo <- FindNeighbors(data.merged.tissue.blood.residents.demo, reduction = "pca", dims = 1:30)
data.merged.tissue.blood.residents.demo <- FindClusters(data.merged.tissue.blood.residents.demo, resolution = 0.6)


saveRDS(data.merged.tissue.blood.residents.demo, "/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/results/data_merged_tissue_blood_residents_DEMO.rds")

# DimPlot(test_2, reduction = "umap", group.by = "final_annotation")

pdf("/home/thatchau/demultiplexing_project/artificial_tissue_blood_experiment/results/UMAP_data_merged_tissue_blood_residents_DEMO.pdf", width = 25, height = 10)

p1.demo <- DimPlot(data.merged.tissue.blood.residents.demo, reduction = "umap", group.by = "orig.ident")
p2.demo <- DimPlot(data.merged.tissue.blood.residents.demo, reduction = "umap", group.by = "final_annotation")

p1.demo + p2.demo

dev.off()