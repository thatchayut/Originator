# This file is used to generate cell type proportion bar plot

set.seed(20211122)

library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

library(speckle)

library(pheatmap)

library(ggplot2)

setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

# Read placenta data with already identified maternal-fetal and blood-tissue immune cell origins
data_MF_TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")

# Remove annotation in which cells are already removed (in this case is `noise`)

data_MF_TB_annotation$annotation = droplevels(data_MF_TB_annotation$annotation, 
                                              exclude = setdiff(levels(data_MF_TB_annotation$annotation),unique(data_MF_TB_annotation$annotation)))

# Remove noise
table(data_MF_TB_annotation$annotation, data_MF_TB_annotation$MF_TB_origin_refWholeblood)

celltype_counts_table_ggplot <- as.data.frame(table(data_MF_TB_annotation$annotation, data_MF_TB_annotation$MF_TB_origin_refWholeblood))

celltype_counts_table_ggplot <- celltype_counts_table_ggplot %>% filter(celltype_counts_table_ggplot$Var1 != "noise")

celltype_counts_table_ggplot$Var1 <- as.character(celltype_counts_table_ggplot$Var1)


# Plot cell type proportion
# Reorder
ref_reorder <- c("Fetal Blood", "Maternal Blood", "Fetal Tissue", "Maternal Tissue")

celltype_counts_table_ggplot_reordered <- celltype_counts_table_ggplot[order(factor(celltype_counts_table_ggplot$Var2, levels = ref_reorder)), ]

# Make Var2 an ordered factor
celltype_counts_table_ggplot_reordered$Var2 <-  factor(x=celltype_counts_table_ggplot_reordered$Var2, 
                                                       levels=ref_reorder)

list_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FFA500", "#800080", 
                 "#FFC0CB", "#00FFFF", "#A52A2A", "#808080", "#008080")

list_ordered_labels_in_plot <- c('Cytotrophoblast', 'Syncytiotrophoblast', 'Extravillous Trophoblast', 'Fibroblast Type 1', 
                                 'Fibroblast Type 2', 'Vascular Endothelial', 'Erythrocyte', 'T-cells',
                                 'Natural Killer cells', 'Macrophage (HB)', 'Monocyte')

ggplot(celltype_counts_table_ggplot_reordered, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = list_colors, 
                    breaks = list_ordered_labels_in_plot) 


pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/celltypeProportion_Hongkong_Placenta_refWholeBlood.pdf", width = 12, height = 10)
ggplot(celltype_counts_table_ggplot_reordered, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values = list_colors, 
                    breaks = list_ordered_labels_in_plot) 
dev.off()

save.image("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/celltype_proportion_placenta_TB_annotation_wholeword.RData")
