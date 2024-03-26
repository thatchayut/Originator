# This file is used to make feature plots of selected genes 

set.seed(20211122)

library(Seurat)

# Read placenta data with already identified maternal-fetal and blood-tissue immune cell origins
data_MF_TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")


### Observe SEPP1 on Macrophage (HB)
data_MF_TB_annotation.MT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Tissue")
data_MF_TB_annotation.FT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Tissue")

data_MF_TB_annotation.macrophage <- subset(data_MF_TB_annotation,
                                           subset = annotation == "Macrophage (HB)")

data_MF_TB_annotation.MT.FT.macrophage <- subset(data_MF_TB_annotation.macrophage,
                                           subset = MF_TB_origin_refWholeblood %in% c("Maternal Tissue", "Fetal Tissue"))

data_MF_TB_annotation.MT.macrophage <- subset(data_MF_TB_annotation.MT, 
                                              subset = annotation == "Macrophage (HB)")

data_MF_TB_annotation.FT.macrophage <- subset(data_MF_TB_annotation.FT, 
                                              subset = annotation == "Macrophage (HB)")

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/FeaturePlot_SEPP1_Hongkong_Placenta_refWholeBlood.pdf", width = 12, height = 10)

DimPlot(data_MF_TB_annotation.MT.FT.macrophage, split.by = "MF_TB_origin_refWholeblood")

FeaturePlot(data_MF_TB_annotation.FT.macrophage, features = c("SEPP1"))
FeaturePlot(data_MF_TB_annotation.MT.macrophage, features = c("SEPP1"))

FeaturePlot(data_MF_TB_annotation.MT.FT.macrophage, features = c("SEPP1"), split.by = "MF_TB_origin_refWholeblood")

dev.off()