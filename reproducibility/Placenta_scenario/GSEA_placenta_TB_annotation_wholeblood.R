# This file is used to perform GSEA comparing common cell types between fetal and maternal tissues

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

setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

# Define a function for getting pathway KEGG
Pathway_KEGG <-function(marker=markers,label=label){
  original_gene_list <- marker$avg_log2FC
  names(original_gene_list) <- row.names(marker)
  ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  df2 = marker[row.names(marker) %in% dedup_ids$SYMBOL,]
  df2$Y = dedup_ids$ENTREZID
  kegg_gene_list <- df2$avg_log2FC
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  GSE <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 3,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 verbose      = FALSE)
  save(GSE, file = paste0(label,"Pathway.RData"))
}

# Read placenta data with already identified maternal-fetal and blood-tissue immune cell origins
data_MF_TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")


# Get common cell types between maternal tissue and fetal tissue
data_MF_TB_annotation.MT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Tissue")
data_MF_TB_annotation.FT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Tissue")

list_celltypes_MT <- unique(data_MF_TB_annotation.MT$annotation)
list_celltypes_FT <- unique(data_MF_TB_annotation.FT$annotation)

ct_list <- intersect(list_celltypes_MT, list_celltypes_FT)
# NOTE: T-cell is excluded as there is only 1 cell present in fetal tissue and is not enough to perform DE
ct_list <- ct_list[ct_list != "T-cells"]

# Get data for each cell prepared for DE analysis

for (ct in ct_list){
  print(ct)
  data_MF_TB_annotation_subset <- subset(x = data_MF_TB_annotation, subset = annotation == ct)
  
  # Compare between fetal and maternal origin 
  Idents(data_MF_TB_annotation_subset) <- "MF_TB_origin_refWholeblood"
  DefaultAssay(data_MF_TB_annotation_subset) <- "RNA"
  markers <- FindMarkers(data_MF_TB_annotation_subset, ident.1 = "Fetal Tissue", ident.2 = "Maternal Tissue", test.use = "MAST", min.pct = 0.1)
  label = paste0(ct, "_Hongkong_FT_MT_")
  # load(file = paste0(label, "DE_genes.RData"))
  try({Pathway_KEGG(marker=markers,label=label)}, silent = FALSE)
  save(markers, file = paste0(label, "DE_genes.RData"))
}

MJ_CT_pathway_table=data.frame()

for (ct in ct_list){
  label = paste0(ct, "_Hongkong_FT_MT_")
  load(file = paste0(label,"Pathway.RData"))
  GSE_tb <- GSE@result
  if(dim(GSE_tb)[1] == 0){
    next
  }
  GSE_tb <- GSE_tb[order(GSE_tb$p.adjust), ]
  if(dim(GSE_tb)[1] >= 30){
    GSE_tb_top30 <- GSE_tb[1:30, ]
  }else{
    GSE_tb_top30 <- GSE_tb
  }
  
  #Temp=data.frame(Pathway=GSE_tb_top30$Description,NES=GSE_tb_top30$NES,Adjusted.Pvalue=GSE_tb_top30$p.adjust,Cell_type=ct)
  Temp=data.frame(id = GSE_tb_top30$ID,Pathway=GSE_tb_top30$Description,NES=GSE_tb_top30$NES,Adjusted.Pvalue=GSE_tb_top30$p.adjust,Cell_type=ct)
  
  MJ_CT_pathway_table=rbind(MJ_CT_pathway_table,Temp)
}

annotation_file <- read.csv("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/kegg_supergroup_anno.csv")

annotation_file_subset <- annotation_file[annotation_file$id %in% unique(MJ_CT_pathway_table$id), ]

subclass_anno <- annotation_file_subset[, c("name", "group", "supergroup")]
subclass_anno <- data.frame(subclass_anno, row.names = 1)

MJ_CT_pathway_table <- merge(MJ_CT_pathway_table, annotation_file_subset, by = "id")

library(tidyr)
library(RColorBrewer)
library(pheatmap)
my.break1 <- seq(-4, 4)
my.color1 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(my.break1))
my.break2 <-my.break1[1:5]
my.color2 <- my.color1[1:5]

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Hongkong_Fetal_Maternal_FetalOriginVMaternalOrigin_KEGG_1.pdf", width = 10, height = 7)
for (i in c("Cellular Processes", "Environmental Information Processing", "Human Diseases", "Genetic Information Processing",
            "Metabolism", "Organismal Systems")){
  temp <- subclass_anno[subclass_anno$supergroup == i, ]
  temp <- temp[, c("group"), drop = F]
  # ptable <- MJ_CT_pathway_table[MJ_CT_pathway_table$Pathway %in% row.names(temp), ]
  ptable <- MJ_CT_pathway_table[MJ_CT_pathway_table$name %in% row.names(temp), ]
  temp <- subclass_anno[subclass_anno$supergroup == i, ]
  
  ptable <- ptable[, c("name", "NES", "Cell_type")]
  ct_dff <- setdiff(ct_list, unique(ptable$Cell_type))
  ptable <- spread(ptable, Cell_type, NES)
  ptable[is.na(ptable)] <- 0
  ptable <- data.frame(ptable, row.names = 1)
  if(length(ct_dff) != 0){
    for(j in ct_dff){
      j <- gsub(" ", ".", j)
      col_ct <- dim(ptable)[2] + 1
      ptable$temp <- rep(0, dim(ptable)[1])
      colnames(ptable)[col_ct] <- j
    }
    ptable <- as.matrix(ptable)
  }else{
    ptable <- as.matrix(ptable)
  }
  if(max(ptable) > 0){
    color_choice = my.color1
    break_choice = my.break1
  }else{
    color_choice = my.color2
    break_choice = my.break2
  }
  pheatmap(ptable, annotation_row = temp, color = color_choice, breaks = break_choice, main = paste0(i))
}
dev.off()

save.image("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/GSEA_placenta_TB_annotation_wholeblood.Rdata")
