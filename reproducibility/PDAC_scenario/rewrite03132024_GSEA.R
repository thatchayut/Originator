# THis file is used to perform Gene Set Enrichment Analysis (GSEA)

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

setwd("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation")

`%!in%` = Negate(`%in%`)

# Read input data
data.TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/PDAC_TB_TuN_annotation_variant_analysis_4_usingWholeblood_03132024.rds")

data.TB_TuN_annotation.cases <- subset(data.TB_TuN_annotation, subset = group == "Cases")


########## Perform GSEA comparing common tissue-resident and blood immune cells in `tumor tissues (cases)`

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

# Get only data immune cells from tumor tissue samples
data.TB_TuN_annotation.immune.cases <- subset(data.TB_TuN_annotation.cases,
                                              subset = compartment == "Immune")

for (ct in list_ct_overlapping_immune_data){
  data.immune_subset <- subset(x = data.TB_TuN_annotation.immune.cases, subset = annotation == ct)
  
  # Compare between this cell type identified as tissue-resident and blood immune cells
  Idents(data.immune_subset) <- "TB_origin_wholeblood"
  DefaultAssay(data.immune_subset) <- "RNA"
  markers <- FindMarkers(data.immune_subset, ident.1 = "Tissue", ident.2 = "Blood", test.use = "wilcox", min.pct = 0.1)
  label = paste0(ct, "_Tissue_Blood_wholebloodref_")
  # load(file = paste0(label, "DE_genes.RData"))
  try({Pathway_KEGG(marker=markers,label=label)}, silent = TRUE)
  save(markers, file = paste0(label, "DE_genes.RData"))
}


MJ_CT_pathway_table=data.frame()

for (ct in list_ct_overlapping_immune_data){
  label = paste0(ct, "_Tissue_Blood_wholebloodref_")
  load(file = paste0(label,"Pathway.RData"))
  GSE_tb <- GSE@result
  if(dim(GSE_tb)[1] == 0){
    next
  }
  GSE_tb <- GSE_tb[order(GSE_tb$p.adjust), ]
  if(dim(GSE_tb)[1] >= 46){
    GSE_tb_top30 <- GSE_tb[1:46, ]
  }else{
    GSE_tb_top30 <- GSE_tb
  }
  
  Temp=data.frame(Pathway=GSE_tb_top30$Description,NES=GSE_tb_top30$NES,Adjusted.Pvalue=GSE_tb_top30$p.adjust,Cell_type=ct)
  MJ_CT_pathway_table=rbind(MJ_CT_pathway_table,Temp)
}


# Get supergroup annotation
annotation_file <- read.csv("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/kegg_supergroup_anno.csv")
annotation_file_subset <- annotation_file[annotation_file$name %in% unique(MJ_CT_pathway_table$Pathway), ]
subclass_anno <- annotation_file_subset[, c("name", "group", "supergroup")]
subclass_anno <- data.frame(subclass_anno, row.names = 1)

# Get heatmap
my.break1 <- seq(-4, 4)
my.color1 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(my.break1))
my.break2 <-my.break1[1:5]
my.color2 <- my.color1[1:5]

unique(subclass_anno$supergroup)

# Without clustering rows
pdf(paste0("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/GSEA_PDAC_TBwholeblood_", "commonImmune", ".SNG.by1000GenomeOverlappingCOSMICvariants_MAFgt0x01Excluded_TumorVsBloodImmune_CasesOnly_KEGG.pdf"), width = 10, height = 7)
for (i in c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
            "Cellular Processes", "Organismal Systems", "Human Diseases")){
  temp <- subclass_anno[subclass_anno$supergroup == i, ]
  temp <- temp[, "group", drop = F]
  ptable <- MJ_CT_pathway_table[MJ_CT_pathway_table$Pathway %in% row.names(temp), ]
  ptable <- ptable[, c("Pathway", "NES", "Cell_type")]
  ct_dff <- setdiff(list_ct_overlapping_immune_data, unique(ptable$Cell_type))
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
  pheatmap(ptable, annotation_row = temp, color = color_choice, breaks = break_choice, main = paste0(i), cluster_rows = FALSE)
}
dev.off()

# With clustering rows
pdf(paste0("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/rewrite011123_results/TB_wholeblood/GSEA_PDAC_TBwholeblood_", "commonImmune", "_TumorVsBloodImmune_CasesOnly_KEGG_clustRows.pdf"), width = 10, height = 7)
for (i in c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
            "Cellular Processes", "Organismal Systems", "Human Diseases")){
  temp <- subclass_anno[subclass_anno$supergroup == i, ]
  temp <- temp[, "group", drop = F]
  ptable <- MJ_CT_pathway_table[MJ_CT_pathway_table$Pathway %in% row.names(temp), ]
  ptable <- ptable[, c("Pathway", "NES", "Cell_type")]
  ct_dff <- setdiff(list_ct_overlapping_immune_data, unique(ptable$Cell_type))
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
  pheatmap(ptable, annotation_row = temp, color = color_choice, breaks = break_choice, main = paste0(i), cluster_rows = TRUE)
}
dev.off()
