# This file is used to perform cell-cell communication (CCC) inference


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

setwd("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation")

`%!in%` = Negate(`%in%`)

# Read input data
data.TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/PDAC_TB_TuN_annotation_variant_analysis_4_usingWholeblood_03132024.rds")

data.TB_annotation.cases <- subset(data.TB_annotation, subset = group == "Cases")

# Prepare the tumor tissue data after removing artifacts
data.TB_annotation.removedArtifacts <- data.TB_annotation

## Remove blood immune cells
data.TB_annotation.removedArtifacts <- subset(data.TB_annotation.removedArtifacts,
                                              subset = TB_origin_wholeblood != "Blood")


# Get only data from cases
data.TB_annotation.removedArtifacts.cases <- subset(data.TB_annotation.removedArtifacts,
                                                    subset = group == "Cases")

########## Perform cellchat in TUMOR TISSUES only comparing BEFORE and AFTER removing artifacts
######## BEFORE: All cell types
######## AFTER: TISSUE immune cells + all other accosiated cells

data.TB_annotation.cases.withArtifacts <- data.TB_annotation.cases
data.TB_annotation.cases.withoutArtifacts <- data.TB_annotation.removedArtifacts.cases

# Set cell type level so that they order the same way in with and without artifacts
ct_level_withArtifacts <- c("Acinar cell", "Beta cell", "Ductal cell", "Endothelial cell", "Epithelial cell",
                            "Fibroblast", "M2 macrophage", "Mast cell", "MDSC", "Monocyte",
                            "NK cell", "Peri-islet Schwann cell", "Plasma cell", "Regulatory T-cell",
                            "T-cell", "B-cell")
ct_level_withoutArtifacts <- c("Acinar cell", "Beta cell", "Ductal cell", "Endothelial cell", "Epithelial cell",
                               "Fibroblast", "M2 macrophage", "Mast cell", "MDSC", "Monocyte",
                               "NK cell", "Peri-islet Schwann cell", "Plasma cell", "Regulatory T-cell",
                               "T-cell")

data.TB_annotation.cases.withArtifacts$annotation <- factor(x = data.TB_annotation.cases.withArtifacts$annotation,
                                                                levels = ct_level_withArtifacts)

data.TB_annotation.cases.withoutArtifacts$annotation <- factor(x = data.TB_annotation.cases.withoutArtifacts$annotation,
                                                                   levels = ct_level_withoutArtifacts)
####### Perform cellchat on tumor tissue WITH artifacts

meta.withArtifacts <- data.TB_annotation.cases.withArtifacts[[]]
cell.use.tt.withArtifacts <- rownames(meta.withArtifacts)

## Remove cell types that were remove earlier from the levels 
meta.withArtifacts$annotation = droplevels(meta.withArtifacts$annotation, 
                                           exclude = setdiff(levels(meta.withArtifacts$annotation),unique(meta.withArtifacts$annotation)))

## Prepare input data for CellChat analysis
data.input.tt.withArtifacts = data.TB_annotation.cases.withArtifacts[, cell.use.tt.withArtifacts]
meta.withArtifacts = meta.withArtifacts[cell.use.tt.withArtifacts, ]

cellchat.tt.withArtifacts <- createCellChat(object = data.input.tt.withArtifacts, meta = meta.withArtifacts, group.by = "annotation")

# Set the ligand-receptor reaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat.tt.withArtifacts@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat.tt.withArtifacts <- subsetData(cellchat.tt.withArtifacts) 

cellchat.tt.withArtifacts <- identifyOverExpressedGenes(cellchat.tt.withArtifacts)
cellchat.tt.withArtifacts <- identifyOverExpressedInteractions(cellchat.tt.withArtifacts)

# Compute the communication probability and infer cellular communication network
cellchat.tt.withArtifacts <- computeCommunProb(cellchat.tt.withArtifacts)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.tt.withArtifacts <- filterCommunication(cellchat.tt.withArtifacts, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.tt.withArtifacts <- computeCommunProbPathway(cellchat.tt.withArtifacts)

# Calculate the aggregated cell-cell communication network
cellchat.tt.withArtifacts <- aggregateNet(cellchat.tt.withArtifacts)

# Visualize the aggregated cell-cell communication network.
## Show the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
groupSize.tt.withArtifacts <- as.numeric(table(cellchat.tt.withArtifacts@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withArtifacts_NumOfInteractions.pdf", width = 12, height = 10)
netVisual_circle(cellchat.tt.withArtifacts@net$count, vertex.weight = groupSize.tt.withArtifacts, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.tt.withArtifacts@net$count, vertex.weight = groupSize.tt.withArtifacts, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.01)

dev.off()

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withArtifacts_InteractionStrength.pdf", width = 12, height = 10)
netVisual_circle(cellchat.tt.withArtifacts@net$weight, vertex.weight = groupSize.tt.withArtifacts, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat.tt.withArtifacts@net$weight, vertex.weight = groupSize.tt.withArtifacts, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.01)
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withArtifacts-individualCelltype.pdf", width = 12, height = 10)
mat.tt.withArtifacts <- cellchat.tt.withArtifacts@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat.tt.withArtifacts)) {
  mat2 <- matrix(0, nrow = nrow(mat.tt.withArtifacts), ncol = ncol(mat.tt.withArtifacts), dimnames = dimnames(mat.tt.withArtifacts))
  mat2[i, ] <- mat.tt.withArtifacts[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.tt.withArtifacts, weight.scale = T, edge.weight.max = max(mat.tt.withArtifacts), title.name = rownames(mat.tt.withArtifacts)[i])
}
dev.off()

# Visualize each signaling pathway
cellchat.tt.withArtifacts@netP$pathways

list_available_pathways <- cellchat.tt.withArtifacts@netP$pathways

for (pathway in list_available_pathways) {
  print(paste0("=> CCC for pathway: ", pathway))
  
  pdf(paste0("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withArtifacts_signalingPathway_", pathway, ".pdf"), width = 12, height = 10)
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "circle")
  
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  
  
  # Chord diagram
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "chord")
  
  # Heatmap
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat.tt.withArtifacts, signaling = pathways.show, color.heatmap = "Reds"))
  
  dev.off()
  
}

saveRDS(cellchat.tt.withArtifacts, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withArtifacts.RDS")

####### Perform cellchat on tumor tissue WITHOUT artifacts
meta.withoutArtifacts <- data.TB_annotation.cases.withoutArtifacts[[]]
cell.use.tt.withoutArtifacts <- rownames(meta.withoutArtifacts)

## Remove cell types that were remove earlier from the levels 
meta.withoutArtifacts$annotation = droplevels(meta.withoutArtifacts$annotation, 
                                              exclude = setdiff(levels(meta.withoutArtifacts$annotation),unique(meta.withoutArtifacts$annotation)))

## Prepare input data for CellChat analysis
data.input.tt.withoutArtifacts = data.TB_annotation.cases.withoutArtifacts[, cell.use.tt.withoutArtifacts]
meta.withoutArtifacts = meta.withoutArtifacts[cell.use.tt.withoutArtifacts, ]

cellchat.tt.withoutArtifacts <- createCellChat(object = data.input.tt.withoutArtifacts, meta = meta.withoutArtifacts, group.by = "annotation")

# Set the ligand-receptor reaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat.tt.withoutArtifacts@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat.tt.withoutArtifacts <- subsetData(cellchat.tt.withoutArtifacts) 

cellchat.tt.withoutArtifacts <- identifyOverExpressedGenes(cellchat.tt.withoutArtifacts)
cellchat.tt.withoutArtifacts <- identifyOverExpressedInteractions(cellchat.tt.withoutArtifacts)

# Compute the communication probability and infer cellular communication network
cellchat.tt.withoutArtifacts <- computeCommunProb(cellchat.tt.withoutArtifacts)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.tt.withoutArtifacts <- filterCommunication(cellchat.tt.withoutArtifacts, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.tt.withoutArtifacts <- computeCommunProbPathway(cellchat.tt.withoutArtifacts)

# Calculate the aggregated cell-cell communication network
cellchat.tt.withoutArtifacts <- aggregateNet(cellchat.tt.withoutArtifacts)

# Visualize the aggregated cell-cell communication network.
## Show the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
groupSize.tt.withoutArtifacts <- as.numeric(table(cellchat.tt.withoutArtifacts@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withoutArtifacts_NumOfInteractions.pdf", width = 12, height = 10)
netVisual_circle(cellchat.tt.withoutArtifacts@net$count, vertex.weight = groupSize.tt.withoutArtifacts, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.tt.withoutArtifacts@net$count, vertex.weight = groupSize.tt.withoutArtifacts, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.01)

dev.off()

pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withoutArtifacts_InteractionStrength.pdf", width = 12, height = 10)
netVisual_circle(cellchat.tt.withoutArtifacts@net$weight, vertex.weight = groupSize.tt.withoutArtifacts, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat.tt.withoutArtifacts@net$weight, vertex.weight = groupSize.tt.withoutArtifacts, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.01)
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks
pdf("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withoutArtifacts-individualCelltype.pdf", width = 12, height = 10)
mat.tt.withoutArtifacts <- cellchat.tt.withoutArtifacts@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat.tt.withoutArtifacts)) {
  mat2 <- matrix(0, nrow = nrow(mat.tt.withoutArtifacts), ncol = ncol(mat.tt.withoutArtifacts), dimnames = dimnames(mat.tt.withoutArtifacts))
  mat2[i, ] <- mat.tt.withoutArtifacts[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.tt.withoutArtifacts, weight.scale = T, edge.weight.max = max(mat.tt.withoutArtifacts), title.name = rownames(mat.tt.withoutArtifacts)[i])
}
dev.off()


# Visualize each signaling pathway
cellchat.tt.withoutArtifacts@netP$pathways

list_available_pathways <- cellchat.tt.withoutArtifacts@netP$pathways

for (pathway in list_available_pathways) {
  print(paste0("=> CCC for pathway: ", pathway))
  
  pdf(paste0("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withoutArtifacts_signalingPathway_", pathway, ".pdf"), width = 12, height = 10)
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "circle")
  
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  
  
  # Chord diagram
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "chord")
  
  # Heatmap
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat.tt.withoutArtifacts, signaling = pathways.show, color.heatmap = "Reds"))
  
  dev.off()
  
}

saveRDS(cellchat.tt.withoutArtifacts, "E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_withoutArtifacts.RDS")

# List overlapping pathway before and after removing artifacts for comparison
list_overlapping_pathways_beforeAfterArtifacts <- intersect(cellchat.tt.withArtifacts@netP$pathways, cellchat.tt.withoutArtifacts@netP$pathways)

for (pathway in list_overlapping_pathways_beforeAfterArtifacts) {
  print(paste0("=> CCC for pathway: ", pathway))
  
  pdf(paste0("E:/Tays_workspace/Demultiplexing_PDAC_new_analysis/new_pipeline/variant_analysis_4/only_TBannotation/results/cellchat/CCC_TumorOnly_compareBeforeAfterRemovingArtifacts_signalingPathway_", pathway, ".pdf"), width = 12, height = 10)
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  p1_withArtifacts <- netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  p1_withoutArtifacts <- netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  p1_withArtifacts
  p1_withoutArtifacts
  
  
  # Circle plot
  par(mfrow=c(1,1))
  p2_withArtifacts <- netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "circle")
  p2_withoutArtifacts <- netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "circle")
  p2_withArtifacts
  p2_withoutArtifacts
  
  
  par(mfrow=c(1,1))
  p3_withArtifacts <- netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  p3_withoutArtifacts <- netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  p3_withArtifacts
  p3_withoutArtifacts
  
  # Chord diagram
  par(mfrow=c(1,1))
  p4_withArtifacts <- netVisual_aggregate(cellchat.tt.withArtifacts, signaling = pathways.show, layout = "chord")
  p4_withoutArtifacts <- netVisual_aggregate(cellchat.tt.withoutArtifacts, signaling = pathways.show, layout = "chord")
  p4_withArtifacts
  p4_withoutArtifacts
  
  # Heatmap
  par(mfrow=c(1,1))
  p5_withArtifacts <- netVisual_heatmap(cellchat.tt.withArtifacts, signaling = pathways.show, color.heatmap = "Reds")
  p5_withoutArtifacts <- netVisual_heatmap(cellchat.tt.withoutArtifacts, signaling = pathways.show, color.heatmap = "Reds")
  
  print(p5_withArtifacts)
  print(p5_withoutArtifacts)
  
  dev.off()
  
}

