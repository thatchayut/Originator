# This file is used to perform DE comparing common cell types between maternal and fetal tissues

library(Seurat)
library(ggVennDiagram)


setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

# Read placenta data with already identified maternal-fetal and blood-tissue immune cell origins
data_MF_TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")

# Get common cell types between maternal tissue and fetal tissue
data_MF_TB_annotation.MT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Tissue")
data_MF_TB_annotation.FT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Tissue")

list_celltypes_MT <- unique(data_MF_TB_annotation.MT$annotation)
list_celltypes_FT <- unique(data_MF_TB_annotation.FT$annotation)

list_commone_celltypes <- intersect(list_celltypes_MT, list_celltypes_FT)

# Get data for each cell prepared for DE analysis
data_MF_TB_annotation.monocyte <- subset(data_MF_TB_annotation, subset = annotation == "Monocyte")
data_MF_TB_annotation.fibroblast1 <- subset(data_MF_TB_annotation, subset = annotation == "Fibroblast Type 1")
data_MF_TB_annotation.fibroblast2 <- subset(data_MF_TB_annotation, subset = annotation == "Fibroblast Type 2")
data_MF_TB_annotation.macrophage <- subset(data_MF_TB_annotation, subset = annotation == "Macrophage (HB)")
data_MF_TB_annotation.vascularEndothelial <- subset(data_MF_TB_annotation, subset = annotation == "Vascular Endothelial")
data_MF_TB_annotation.tcell <- subset(data_MF_TB_annotation, subset = annotation == "T-cells")
data_MF_TB_annotation.nk <- subset(data_MF_TB_annotation, subset = annotation == "Natural Killer cells")

# Perform DE
Idents(data_MF_TB_annotation.monocyte) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.fibroblast1) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.fibroblast2) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.macrophage) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.vascularEndothelial) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.tcell) <- "MF_TB_origin_refWholeblood"
Idents(data_MF_TB_annotation.nk) <- "MF_TB_origin_refWholeblood"

DefaultAssay(data_MF_TB_annotation.monocyte) <- "RNA"
DefaultAssay(data_MF_TB_annotation.fibroblast1) <- "RNA"
DefaultAssay(data_MF_TB_annotation.fibroblast2) <- "RNA"
DefaultAssay(data_MF_TB_annotation.macrophage) <- "RNA"
DefaultAssay(data_MF_TB_annotation.vascularEndothelial) <- "RNA"
DefaultAssay(data_MF_TB_annotation.tcell) <- "RNA"
DefaultAssay(data_MF_TB_annotation.nk) <- "RNA"

DE.monocyte.FT.MT <- FindMarkers(data_MF_TB_annotation.monocyte, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
DE.fibroblast1.FT.MT <- FindMarkers(data_MF_TB_annotation.fibroblast1, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
DE.fibroblast2.FT.MT <- FindMarkers(data_MF_TB_annotation.fibroblast2, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
DE.macrophage.FT.MT <- FindMarkers(data_MF_TB_annotation.macrophage, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
DE.vascularEndothelial.FT.MT <- FindMarkers(data_MF_TB_annotation.vascularEndothelial, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
DE.nk.FT.MT <- FindMarkers(data_MF_TB_annotation.nk, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", only.pos = FALSE, test.use = "wilcox")
# NOTE: T-cell is excluded as there is only 1 cell present in fetal tissue and is not enough to perform DE

# Sort by adj. p-value (ascending) and log2fc (descending)
DE.monocyte.FT.MT <- DE.monocyte.FT.MT[order(-DE.monocyte.FT.MT$avg_log2FC,
                                             DE.monocyte.FT.MT$p_val_adj),]

DE.fibroblast1.FT.MT <- DE.fibroblast1.FT.MT[order(-DE.fibroblast1.FT.MT$avg_log2FC,
                                                   DE.fibroblast1.FT.MT$p_val_adj),]

DE.fibroblast2.FT.MT <- DE.fibroblast2.FT.MT[order(-DE.fibroblast2.FT.MT$avg_log2FC,
                                                   DE.fibroblast2.FT.MT$p_val_adj),]

DE.macrophage.FT.MT <- DE.macrophage.FT.MT[order(-DE.macrophage.FT.MT$avg_log2FC,
                                                 DE.macrophage.FT.MT$p_val_adj),]

DE.vascularEndothelial.FT.MT <- DE.vascularEndothelial.FT.MT[order(-DE.vascularEndothelial.FT.MT$avg_log2FC,
                                                                   DE.vascularEndothelial.FT.MT$p_val_adj),]

DE.nk.FT.MT <- DE.nk.FT.MT[order(-DE.nk.FT.MT$avg_log2FC,
                                 DE.nk.FT.MT$p_val_adj),]

# Filter to get DEG with p-value < significance
significance <- 0.01

DE.monocyte.FT.MT_significance <- DE.monocyte.FT.MT %>% filter(p_val_adj < significance)
DE.fibroblast1.FT.MT_significance <- DE.fibroblast1.FT.MT %>% filter(p_val_adj < significance)
DE.fibroblast2.FT.MT_significance <- DE.fibroblast2.FT.MT %>% filter(p_val_adj < significance)
DE.macrophage.FT.MT_significance <- DE.macrophage.FT.MT %>% filter(p_val_adj < significance)
DE.vascularEndothelial.FT.MT_significance <- DE.vascularEndothelial.FT.MT %>% filter(p_val_adj < significance)
DE.nk.FT.MT_significance <- DE.nk.FT.MT %>% filter(p_val_adj < significance)

# Save DEGs
saveRDS(DE.monocyte.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_monocyte_FT-MT_significance.rds")
write.csv(DE.monocyte.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_monocyte_FT-MT_significance.csv")

saveRDS(DE.fibroblast1.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast1_FT-MT_significance.rds")
write.csv(DE.fibroblast1.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast1_FT-MT_significance.csv")

saveRDS(DE.fibroblast2.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast2_FT-MT_significance.rds")
write.csv(DE.fibroblast2.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast2_FT-MT_significance.csv")

saveRDS(DE.macrophage.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_macrophage_FT-MT_significance.rds")
write.csv(DE.macrophage.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_macrophage_FT-MT_significance.csv")

saveRDS(DE.vascularEndothelial.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_vascularEndothelial_FT-MT_significance.rds")
write.csv(DE.vascularEndothelial.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_vascularEndothelial_FT-MT_significance.csv")

saveRDS(DE.nk.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_NKcell_FT-MT_significance.rds")
write.csv(DE.nk.FT.MT_significance, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_NKcell_FT-MT_significance.csv")


# Prepare data frame for Venn diagram
DE.monocyte.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_monocyte_FT-MT_significance.rds")
DE.fibroblast1.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast1_FT-MT_significance.rds")
DE.fibroblast2.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast2_FT-MT_significance.rds")
DE.macrophage.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_macrophage_FT-MT_significance.rds")
DE.vascularEndothelial.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_vascularEndothelial_FT-MT_significance.rds")
DE.nk.FT.MT_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_NKcell_FT-MT_significance.rds")


df_DE.monocyte.FT.MT_significance <- data.frame(cell_type = rep("Monocyte",
                                                               nrow(DE.monocyte.FT.MT_significance)),
                                                       DEGs = row.names(DE.monocyte.FT.MT_significance))

df_DE.fibroblast1.FT.MT_significance <- data.frame(cell_type = rep("Fibroblast type 1",
                                                               nrow(DE.fibroblast1.FT.MT_significance)),
                                               DEGs = row.names(DE.fibroblast1.FT.MT_significance))

df_DE.fibroblast2.FT.MT_significance <- data.frame(cell_type = rep("Fibroblast type 2",
                                                                   nrow(DE.fibroblast2.FT.MT_significance)),
                                                   DEGs = row.names(DE.fibroblast2.FT.MT_significance))

df_DE.macrophage.FT.MT_significance <- data.frame(cell_type = rep("Macrophage (HB)",
                                                                   nrow(DE.macrophage.FT.MT_significance)),
                                                   DEGs = row.names(DE.macrophage.FT.MT_significance))

df_DE.vascularEndothelial.FT.MT_significance <- data.frame(cell_type = rep("Vascular endothelial",
                                                                          nrow(DE.vascularEndothelial.FT.MT_significance)),
                                                  DEGs = row.names(DE.vascularEndothelial.FT.MT_significance))

df_DE.nk.FT.MT_significance <- data.frame(cell_type = rep("NK cell",
                                                           nrow(DE.nk.FT.MT_significance)),
                                                           DEGs = row.names(DE.nk.FT.MT_significance))

# NOTE: Disregard B-cell as no tissue-resident B-cell identified

list_common_ct_de_markers_all_log2FC_significance <- list(monocyte = df_DE.monocyte.FT.MT_significance$DEGs,
                                                         fibroblast1 = df_DE.fibroblast1.FT.MT_significance$DEGs,
                                                         fibroblast2 = df_DE.fibroblast2.FT.MT_significance$DEGs,
                                                         macrophage = df_DE.macrophage.FT.MT_significance$DEGs,
                                                         vascular_endothelial = df_DE.vascularEndothelial.FT.MT_significance$DEGs,
                                                         NK = df_DE.nk.FT.MT_significance$DEGs)

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/plaenta_TBwholeblood_VennDiagram_DEGs_commonCelltype_FTvsMT.pdf", width = 12, height = 10)
ggVennDiagram(list_common_ct_de_markers_all_log2FC_significance, category.names = c("Monocyte",
                                                                                     "Fibroblast type 1",
                                                                                     "Fibroblast type 2",
                                                                                     "Macrophage",
                                                                                     "Vascular endothelial cell",
                                                                                    "NK cell"))
ggVennDiagram(list_common_ct_de_markers_all_log2FC_significance, category.names = c("Monocyte",
                                                                                    "Fibroblast type 1",
                                                                                    "Fibroblast type 2",
                                                                                    "Macrophage",
                                                                                    "Vascular endothelial cell",
                                                                                    "NK cell"),
              label = "none")

dev.off()


