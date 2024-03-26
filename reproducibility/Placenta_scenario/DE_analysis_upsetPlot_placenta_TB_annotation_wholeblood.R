# This file is used to make an upset plot (alternative of Venn diagram) of DEGs comparing
#   fetal and maternal tissues

set.seed(20211122)
library(ggplot2);
library(data.table);
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyr)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library("RColorBrewer")
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

library(clustree)

library(tidyr)
library(RColorBrewer)
library(pheatmap)

library(CellChat)

library(ggVennDiagram)

library(SeuratDisk)

library(caret)

library(UpSetR)

setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

`%!in%` = Negate(`%in%`)

# Prepare DEGs
DE.immune.cases.fibroblast1.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast1_FT-MT_significance.rds")
DE.immune.cases.fibroblast2.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_fibroblast2_FT-MT_significance.rds")
DE.immune.cases.macrophage.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_macrophage_FT-MT_significance.rds")
DE.immune.cases.vascularEndothelial.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_vascularEndothelial_FT-MT_significance.rds")
DE.immune.cases.NK.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_NKcell_FT-MT_significance.rds")
DE.immune.cases.monocyte.TB_significance <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_immune_monocyte_FT-MT_significance.rds")

# Prepare info of intersected DEGs
list_ct <- c("fibroblast1", "fibroblast2", "macrophage", "vascular endothelial", "NK cell", "monocyte")


## Intersection between DEGs in 2 cell types
combn(list_ct, 2)

upset_2ct_fibroblast1_fibroblast2 <- intersect(row.names(DE.immune.cases.fibroblast1.TB_significance),
                                               row.names(DE.immune.cases.fibroblast2.TB_significance))

upset_2ct_fibroblast1_macrophage <- intersect(row.names(DE.immune.cases.fibroblast1.TB_significance),
                                              row.names(DE.immune.cases.macrophage.TB_significance))

upset_2ct_fibroblast1_vascularEndothelial <- intersect(row.names(DE.immune.cases.fibroblast1.TB_significance),
                                                       row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_2ct_fibroblast1_NK <- intersect(row.names(DE.immune.cases.fibroblast1.TB_significance),
                                      row.names(DE.immune.cases.NK.TB_significance))

upset_2ct_fibroblast1_monocyte <- intersect(row.names(DE.immune.cases.fibroblast1.TB_significance),
                                            row.names(DE.immune.cases.monocyte.TB_significance))

upset_2ct_fibroblast2_macrophage <- intersect(row.names(DE.immune.cases.fibroblast2.TB_significance),
                                              row.names(DE.immune.cases.macrophage.TB_significance))

upset_2ct_fibroblast2_vascularEndothelial <- intersect(row.names(DE.immune.cases.fibroblast2.TB_significance),
                                                       row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_2ct_fibroblast2_NK <- intersect(row.names(DE.immune.cases.fibroblast2.TB_significance),
                                      row.names(DE.immune.cases.NK.TB_significance))

upset_2ct_fibroblast2_monocyte <- intersect(row.names(DE.immune.cases.fibroblast2.TB_significance),
                                            row.names(DE.immune.cases.monocyte.TB_significance))

upset_2ct_macrophage_vascularEndothelial <- intersect(row.names(DE.immune.cases.macrophage.TB_significance),
                                                      row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_2ct_macrophage_NK <- intersect(row.names(DE.immune.cases.macrophage.TB_significance),
                                     row.names(DE.immune.cases.NK.TB_significance))

upset_2ct_macrophage_monocyte <- intersect(row.names(DE.immune.cases.macrophage.TB_significance),
                                           row.names(DE.immune.cases.monocyte.TB_significance))

upset_2ct_vascularEndothelial_NK <- intersect(row.names(DE.immune.cases.vascularEndothelial.TB_significance),
                                              row.names(DE.immune.cases.NK.TB_significance))

upset_2ct_vascularEndothelial_monocyte <- intersect(row.names(DE.immune.cases.vascularEndothelial.TB_significance),
                                                    row.names(DE.immune.cases.monocyte.TB_significance))

upset_2ct_NK_monocyte <- intersect(row.names(DE.immune.cases.NK.TB_significance),
                                   row.names(DE.immune.cases.monocyte.TB_significance))

## Intersection between DEGs in 3 cell types
combn(list_ct, 3)

upset_3ct_fibroblast1_fibroblast2_macrophage <- intersect(upset_2ct_fibroblast1_fibroblast2, 
                                                          row.names(DE.immune.cases.macrophage.TB_significance))

upset_3ct_fibroblast1_fibroblast2_vascularEndothelial <- intersect(upset_2ct_fibroblast1_fibroblast2, 
                                                                  row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_3ct_fibroblast1_fibroblast2_NK <- intersect(upset_2ct_fibroblast1_fibroblast2, 
                                                  row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_fibroblast1_fibroblast2_monocyte <- intersect(upset_2ct_fibroblast1_fibroblast2, 
                                                        row.names(upset_2ct_fibroblast1_monocyte))

upset_3ct_fibroblast1_macrophage_vascularEndothelial <- intersect(upset_2ct_fibroblast1_macrophage, 
                                                                  row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_3ct_fibroblast1_macrophage_NK <- intersect(upset_2ct_fibroblast1_macrophage, 
                                                 row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_fibroblast1_macrophage_monocyte <- intersect(upset_2ct_fibroblast1_macrophage, 
                                                       row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_fibroblast1_vascularEndothelial_NK <- intersect(upset_2ct_fibroblast1_vascularEndothelial, 
                                                          row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_fibroblast1_vascularEndothelial_monocyte <- intersect(upset_2ct_fibroblast1_vascularEndothelial, 
                                                                row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_fibroblast1_NK_monocyte <- intersect(upset_2ct_fibroblast1_NK, 
                                               row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_fibroblast2_macrophage_vascularEndothelial <- intersect(upset_2ct_fibroblast2_macrophage, 
                                                                  row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_3ct_fibroblast2_macrophage_NK <- intersect(upset_2ct_fibroblast2_macrophage, 
                                                 row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_fibroblast2_macrophage_monocyte <- intersect(upset_2ct_fibroblast2_macrophage, 
                                                       row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_fibroblast2_vascularEndothelial_NK <- intersect(upset_2ct_fibroblast2_vascularEndothelial, 
                                                          row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_fibroblast2_vascularEndothelial_monocyte <- intersect(upset_2ct_fibroblast2_vascularEndothelial, 
                                                                row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_fibroblast2_NK_monocyte <- intersect(upset_2ct_fibroblast2_NK, 
                                               row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_macrophage_vascularEndothelial_NK <- intersect(upset_2ct_macrophage_vascularEndothelial, 
                                                         row.names(DE.immune.cases.NK.TB_significance))

upset_3ct_macrophage_vascularEndothelial_monocyte <- intersect(upset_2ct_macrophage_vascularEndothelial, 
                                                               row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_macrophage_NK_monocyte <- intersect(upset_2ct_macrophage_NK, 
                                              row.names(DE.immune.cases.monocyte.TB_significance))

upset_3ct_vascularEndothelial_NK_monocyte <- intersect(upset_2ct_vascularEndothelial_NK, 
                                                       row.names(DE.immune.cases.monocyte.TB_significance))

## Intersection between DEGs in 4 cell types
combn(list_ct, 4)

upset_4ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial <- intersect(upset_3ct_fibroblast1_fibroblast2_macrophage, 
                                                                              row.names(DE.immune.cases.vascularEndothelial.TB_significance))

upset_4ct_fibroblast1_fibroblast2_macrophage_NK<- intersect(upset_3ct_fibroblast1_fibroblast2_macrophage, 
                                                            row.names(DE.immune.cases.NK.TB_significance))

upset_4ct_fibroblast1_fibroblast2_macrophage_monocyte<- intersect(upset_3ct_fibroblast1_fibroblast2_macrophage, 
                                                                  row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast1_fibroblast2_vascularEndothelial_NK <- intersect(upset_3ct_fibroblast1_fibroblast2_vascularEndothelial, 
                                                                      row.names(DE.immune.cases.NK.TB_significance))

upset_4ct_fibroblast1_fibroblast2_vascularEndothelial_monocyte <- intersect(upset_3ct_fibroblast1_fibroblast2_vascularEndothelial, 
                                                                            row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast1_fibroblast2_NK_monocyte <- intersect(upset_3ct_fibroblast1_fibroblast2_NK, 
                                                           row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast1_macrophage_vascularEndothelial_NK <- intersect(upset_3ct_fibroblast1_macrophage_vascularEndothelial, 
                                                                     row.names(DE.immune.cases.NK.TB_significance))

upset_4ct_fibroblast1_macrophage_vascularEndothelial_monocyte <- intersect(upset_3ct_fibroblast1_macrophage_vascularEndothelial, 
                                                                           row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast1_macrophage_NK_monocyte <- intersect(upset_3ct_fibroblast1_macrophage_NK, 
                                                          row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast1_vascularEndothelial_NK_monocyte <- intersect(upset_3ct_fibroblast1_vascularEndothelial_NK, 
                                                                   row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast2_macrophage_vascularEndothelial_NK <- intersect(upset_3ct_fibroblast2_macrophage_vascularEndothelial, 
                                                                     row.names(DE.immune.cases.NK.TB_significance))

upset_4ct_fibroblast2_macrophage_vascularEndothelial_monocyte <- intersect(upset_3ct_fibroblast2_macrophage_vascularEndothelial, 
                                                                           row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast2_macrophage_NK_monocyte <- intersect(upset_3ct_fibroblast2_macrophage_NK, 
                                                          row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_fibroblast2_vascularEndothelial_NK_monocyte <- intersect(upset_3ct_fibroblast2_vascularEndothelial_NK, 
                                                                   row.names(DE.immune.cases.monocyte.TB_significance))

upset_4ct_macrophage_vascularEndothelial_NK_monocyte <- intersect(upset_3ct_macrophage_vascularEndothelial_NK, 
                                                                  row.names(DE.immune.cases.monocyte.TB_significance))

## Intersection between DEGs in 4 cell types
combn(list_ct, 5)

upset_5ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_NK <- intersect(upset_4ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial, 
                                                                              row.names(DE.immune.cases.NK.TB_significance))

upset_5ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_monocyte <- intersect(upset_4ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial, 
                                                                                      row.names(DE.immune.cases.monocyte.TB_significance))

upset_5ct_fibroblast1_fibroblast2_macrophage_NK_monocyte <- intersect(upset_4ct_fibroblast1_fibroblast2_macrophage_NK, 
                                                                      row.names(DE.immune.cases.monocyte.TB_significance))

upset_5ct_fibroblast1_fibroblast2_vascularEndothelial_NK_monocyte <- intersect(upset_4ct_fibroblast1_fibroblast2_vascularEndothelial_NK, 
                                                                      row.names(DE.immune.cases.monocyte.TB_significance))

upset_5ct_fibroblast1_macrophage_vascularEndothelial_NK_monocyte <- intersect(upset_4ct_fibroblast1_macrophage_vascularEndothelial_NK, 
                                                                               row.names(DE.immune.cases.monocyte.TB_significance))

upset_5ct_fibroblast2_macrophage_vascularEndothelial_NK_monocyte <- intersect(upset_4ct_fibroblast2_macrophage_vascularEndothelial_NK, 
                                                                              row.names(DE.immune.cases.monocyte.TB_significance))

## Intersection between DEGs in 4 cell types
combn(list_ct, 6)

upset_6ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_NK_monocyte <- intersect(upset_5ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_NK, 
                                                                                          row.names(DE.immune.cases.monocyte.TB_significance))

# Prepare data for UpsetR plot
upsetR_input_DE <- c(fibroblast1 = nrow(DE.immune.cases.fibroblast1.TB_significance),
                     fibroblast2 = nrow(DE.immune.cases.fibroblast2.TB_significance),
                     macrophage = nrow(DE.immune.cases.macrophage.TB_significance),
                     vascular.endothelial = nrow(DE.immune.cases.vascularEndothelial.TB_significance),
                     NK = nrow(DE.immune.cases.NK.TB_significance),
                     monocyte = nrow(DE.immune.cases.monocyte.TB_significance),
                     
                     ## Intersection between DEGs in 2 cell types
                     "fibroblast1&fibroblast2" = length(upset_2ct_fibroblast1_fibroblast2),
                     "fibroblast1&macrophage" = length(upset_2ct_fibroblast1_macrophage),
                     "fibroblast1&vascular.endothelial" = length(upset_2ct_fibroblast1_vascularEndothelial),
                     "fibroblast1&NK" = length(upset_2ct_fibroblast1_NK),
                     "fibroblast1&monocyte" = length(upset_2ct_fibroblast1_monocyte),
                     "fibroblast2&macrophage" = length(upset_2ct_fibroblast2_macrophage),
                     "fibroblast2&vascular.endothelial" = length(upset_2ct_fibroblast2_vascularEndothelial),
                     "fibroblast2&NK" = length(upset_2ct_fibroblast2_NK),
                     "fibroblast2&monocyte" = length(upset_2ct_fibroblast2_monocyte),
                     "macrophage&vascular.endothelial" = length(upset_2ct_macrophage_vascularEndothelial),
                     "macrophage&NK" = length(upset_2ct_macrophage_NK),
                     "macrophage&monocyte" = length(upset_2ct_macrophage_monocyte),
                     "vascular.endothelial&NK" = length(upset_2ct_vascularEndothelial_NK),
                     "vascular.endothelial&monocyte" = length(upset_2ct_vascularEndothelial_monocyte),
                     "NK&monocyte" = length(upset_2ct_NK_monocyte),
                     
                     ## Intersection between DEGs in 3 cell types
                     "fibroblast1&fibroblast2&macrophage" = length(upset_3ct_fibroblast1_fibroblast2_macrophage),
                     "fibroblast1&fibroblast2&vascular.endothelial" = length(upset_3ct_fibroblast1_fibroblast2_vascularEndothelial),
                     "fibroblast1&fibroblast2&NK" = length(upset_3ct_fibroblast1_fibroblast2_NK),
                     "fibroblast1&fibroblast2&monocyte" = length(upset_3ct_fibroblast1_fibroblast2_monocyte),
                     "fibroblast1&macrophage&vascular.endothelial" = length(upset_3ct_fibroblast1_macrophage_vascularEndothelial),
                     "fibroblast1&macrophage&NK" = length(upset_3ct_fibroblast1_macrophage_NK),
                     "fibroblast1&macrophage&monocyte" = length(upset_3ct_fibroblast1_macrophage_monocyte),
                     "fibroblast1&vascular.endothelial&NK" = length(upset_3ct_fibroblast1_vascularEndothelial_NK),
                     "fibroblast1&vascular.endothelial&monocyte" = length(upset_3ct_fibroblast1_vascularEndothelial_monocyte),
                     "fibroblast1&NK&monocyte" = length(upset_3ct_fibroblast1_NK_monocyte),
                     "fibroblast2&macrophage&vascular.endothelial" = length(upset_3ct_fibroblast2_macrophage_vascularEndothelial),
                     "fibroblast2&macrophage&NK" = length(upset_3ct_fibroblast2_macrophage_NK),
                     "fibroblast2&macrophage&monocyte" = length(upset_3ct_fibroblast2_macrophage_monocyte),
                     "fibroblast2&vascular.endothelial&NK" = length(upset_3ct_fibroblast2_vascularEndothelial_NK),
                     "fibroblast2&vascular.endothelial&monocyte" = length(upset_3ct_fibroblast2_vascularEndothelial_monocyte),
                     "fibroblast2&NK&monocyte" = length(upset_3ct_fibroblast2_NK_monocyte),
                     "macrophage&vascular.endothelial&NK" = length(upset_3ct_macrophage_vascularEndothelial_NK),
                     "macrophage&vascular.endothelial&monocyte" = length(upset_3ct_macrophage_vascularEndothelial_monocyte),
                     "macrophage&NK&monocyte" = length(upset_3ct_macrophage_NK_monocyte),
                     "vascular.endothelial&NK&monocyte" = length(upset_3ct_vascularEndothelial_NK_monocyte),
                     
                     ## Intersection between DEGs in 4 cell types
                     "fibroblast1&fibroblast2&macrophage&vascular.endothelial" = length(upset_4ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial),
                     "fibroblast1&fibroblast2&macrophage&NK" = length(upset_4ct_fibroblast1_fibroblast2_macrophage_NK),
                     "fibroblast1&fibroblast2&macrophage&monocyte" = length(upset_4ct_fibroblast1_fibroblast2_macrophage_monocyte),
                     "fibroblast1&fibroblast2&vascular.endothelial&NK" = length(upset_4ct_fibroblast1_fibroblast2_vascularEndothelial_NK),
                     "fibroblast1&fibroblast2&vascular.endothelial&monocyte" = length(upset_4ct_fibroblast1_fibroblast2_vascularEndothelial_monocyte),
                     "fibroblast1&fibroblast2&NK&monocyte" = length(upset_4ct_fibroblast1_fibroblast2_NK_monocyte),
                     "fibroblast1&macrophage&vascular.endothelial&NK" = length(upset_4ct_fibroblast1_macrophage_vascularEndothelial_NK),
                     "fibroblast1&macrophage&vascular.endothelial&monocyte" = length(upset_4ct_fibroblast1_macrophage_vascularEndothelial_monocyte),
                     "fibroblast1&macrophage&NK&monocyte" = length(upset_4ct_fibroblast1_macrophage_NK_monocyte),
                     "fibroblast1&macrophage&vascular.endothelial&NK" = length(upset_4ct_fibroblast1_macrophage_vascularEndothelial_NK),
                     "fibroblast1&macrophage&vascular.endothelial&monocyte" = length(upset_4ct_fibroblast1_macrophage_vascularEndothelial_monocyte),
                     "fibroblast1&macrophage&NK&monocyte" = length(upset_4ct_fibroblast1_macrophage_NK_monocyte),
                     "fibroblast2&macrophage&vascular.endothelial&NK" = length(upset_4ct_fibroblast2_macrophage_vascularEndothelial_NK),
                     "fibroblast2&macrophage&vascular.endothelial&monocyte" = length(upset_4ct_fibroblast2_macrophage_vascularEndothelial_monocyte),
                     "fibroblast2&macrophage&NK&monocyte" = length(upset_4ct_fibroblast2_macrophage_NK_monocyte),
                     "fibroblast2&vascular.endothelial&NK&monocyte" = length(upset_4ct_fibroblast2_vascularEndothelial_NK_monocyte),
                     "macrophage&vascular.endothelial&NK&monocyte" = length(upset_4ct_macrophage_vascularEndothelial_NK_monocyte),
                     
                     ## Intersection between DEGs in 5 cell types
                     "fibroblast1&fibroblast2&macrophage&vascular.endothelial&NK" = length(upset_5ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_NK),
                     "fibroblast1&fibroblast2&macrophage&vascular.endothelial&monocyte" = length(upset_5ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_monocyte),
                     "fibroblast1&fibroblast2&macrophage&NK&monocyte" = length(upset_5ct_fibroblast1_fibroblast2_macrophage_NK_monocyte),
                     "fibroblast1&fibroblast2&vascular.endothelial&NK&monocyte" = length(upset_5ct_fibroblast1_fibroblast2_vascularEndothelial_NK_monocyte),
                     "fibroblast1&macrophage&vascular.endothelial&NK&monocyte" = length(upset_5ct_fibroblast1_macrophage_vascularEndothelial_NK_monocyte),
                     "fibroblast2&macrophage&vascular.endothelial&NK&monocyte" = length(upset_5ct_fibroblast2_macrophage_vascularEndothelial_NK_monocyte),
                     
                     ## Intersection between DEGs in 6 cell types
                     "fibroblast1&fibroblast2&macrophage&vascular.endothelial&NK&monocyte" = length(upset_6ct_fibroblast1_fibroblast2_macrophage_vascularEndothelial_NK_monocyte)
                    )
                     
pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TBwholeblood_UpsetPlot_DEGs_commonCelltType_FTvsMT.pdf", width = 12, height = 10)
upset(fromExpression(upsetR_input_DE), order.by = "freq", decreasing = TRUE, nsets = 6)
dev.off()

save.image("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/DE_analysis_upsetPlot_placenta_TB_annotation_wholeblood.RData")
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     

