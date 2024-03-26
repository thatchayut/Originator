# Perform cell-cell communication inference on fetal tissue and maternal tissue

set.seed(20211122)

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

setwd("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/")

# Read placenta data with already identified maternal-fetal and blood-tissue immune cell origins
data_MF_TB_annotation <- readRDS("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/Placenta_TB_MF_annotation_usingWholeblood_022424.rds")

# Get common cell types between maternal tissue and fetal tissue
data_MF_TB_annotation.MT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Maternal Tissue")
data_MF_TB_annotation.FT <- subset(data_MF_TB_annotation, subset = MF_TB_origin_refWholeblood == "Fetal Tissue")

list_celltypes_MT <- unique(data_MF_TB_annotation.MT$annotation)
list_celltypes_FT <- unique(data_MF_TB_annotation.FT$annotation)

list_commone_celltypes <- intersect(list_celltypes_MT, list_celltypes_FT)

# Filter to get only common cell types
data_MF_TB_annotation.MT.commonCelltype <- subset(data_MF_TB_annotation.MT,
                                                  subset = annotation %in%  list_commone_celltypes)
data_MF_TB_annotation.FT.commonCelltype <- subset(data_MF_TB_annotation.FT,
                                                  subset = annotation %in%  list_commone_celltypes)

DefaultAssay(data_MF_TB_annotation.MT.commonCelltype) <- "RNA"
DefaultAssay(data_MF_TB_annotation.FT.commonCelltype) <- "RNA"

######### Maternal tissue
meta.MT <- data_MF_TB_annotation.MT.commonCelltype[[]]
cell.use.MT <- rownames(meta.MT)

## Remove cell types that were remove earlier from the levels 
meta.MT$annotation = droplevels(meta.MT$annotation, 
                                   exclude = setdiff(levels(meta.MT$annotation),unique(meta.MT$annotation)))

## Prepare input data for CellChat analysis
data_MF_TB_annotation.MT.commonCelltype.cellchat = data_MF_TB_annotation.MT.commonCelltype[, cell.use.MT]
meta.MT = meta.MT[cell.use.MT, ]

cellchat.MT <- createCellChat(object = data_MF_TB_annotation.MT.commonCelltype.cellchat, meta = meta.MT, group.by = "annotation")

# Set the ligand-receptor reaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat.MT@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat.MT <- subsetData(cellchat.MT) 

cellchat.MT <- identifyOverExpressedGenes(cellchat.MT)
cellchat.MT <- identifyOverExpressedInteractions(cellchat.MT)

# Compute the communication probability and infer cellular communication network
cellchat.MT <- computeCommunProb(cellchat.MT)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.MT <- filterCommunication(cellchat.MT, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.MT <- computeCommunProbPathway(cellchat.MT)

# Calculate the aggregated cell-cell communication network
cellchat.MT <- aggregateNet(cellchat.MT)

# Visualize the aggregated cell-cell communication network.
## Show the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
groupSize.MT <- as.numeric(table(cellchat.MT@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_MaternalTissue_NumOfInteractions.pdf", width = 12, height = 10)
netVisual_circle(cellchat.MT@net$count, vertex.weight = groupSize.MT, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.MT@net$count, vertex.weight = groupSize.MT, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.01)

dev.off()

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_MaternalTissue_withArtifacts_InteractionStrength.pdf", width = 12, height = 10)
netVisual_circle(cellchat.MT@net$weight, vertex.weight = groupSize.MT, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat.MT@net$weight, vertex.weight = groupSize.MT, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.01)
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks
pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_MaternalTissue-individualCelltype.pdf", width = 12, height = 10)
mat.MT <- cellchat.MT@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat.MT)) {
  mat2 <- matrix(0, nrow = nrow(mat.MT), ncol = ncol(mat.MT), dimnames = dimnames(mat.MT))
  mat2[i, ] <- mat.MT[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.MT, weight.scale = T, edge.weight.max = max(mat.MT), title.name = rownames(mat.MT)[i])
}
dev.off()

# Visualize each signaling pathway
cellchat.MT@netP$pathways

list_available_pathways <- cellchat.MT@netP$pathways

for (pathway in list_available_pathways) {
  print(paste0("=> CCC for pathway: ", pathway))
  
  pdf(paste0("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_MaternalTissue_signalingPathway_", pathway, ".pdf"), width = 12, height = 10)
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat.MT, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.MT, signaling = pathways.show, layout = "circle")
  
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.MT, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  
  
  # Chord diagram
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.MT, signaling = pathways.show, layout = "chord")
  
  # Heatmap
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat.MT, signaling = pathways.show, color.heatmap = "Reds"))
  
  dev.off()
  
}

saveRDS(cellchat.MT, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_MaternalTissue.RDS")

######### Fetal tissue
meta.FT <- data_MF_TB_annotation.FT.commonCelltype[[]]
cell.use.FT <- rownames(meta.FT)

## Remove cell types that were remove earlier from the levels 
meta.FT$annotation = droplevels(meta.FT$annotation, 
                                exclude = setdiff(levels(meta.FT$annotation),unique(meta.FT$annotation)))

## Prepare input data for CellChat analysis
data_MF_TB_annotation.FT.commonCelltype.cellchat = data_MF_TB_annotation.FT.commonCelltype[, cell.use.FT]
meta.FT = meta.FT[cell.use.FT, ]

cellchat.FT <- createCellChat(object = data_MF_TB_annotation.FT.commonCelltype.cellchat, meta = meta.FT, group.by = "annotation")

# Set the ligand-receptor reaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB 

# set the used database in the object
cellchat.FT@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat.FT <- subsetData(cellchat.FT) 

cellchat.FT <- identifyOverExpressedGenes(cellchat.FT)
cellchat.FT <- identifyOverExpressedInteractions(cellchat.FT)

# Compute the communication probability and infer cellular communication network
cellchat.FT <- computeCommunProb(cellchat.FT)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.FT <- filterCommunication(cellchat.FT, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat.FT <- computeCommunProbPathway(cellchat.FT)

# Calculate the aggregated cell-cell communication network
cellchat.FT <- aggregateNet(cellchat.FT)

# Visualize the aggregated cell-cell communication network.
## Show the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot
groupSize.FT <- as.numeric(table(cellchat.FT@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_FetalTissue_NumOfInteractions.pdf", width = 12, height = 10)
netVisual_circle(cellchat.FT@net$count, vertex.weight = groupSize.FT, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.FT@net$count, vertex.weight = groupSize.FT, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.01)

dev.off()

pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_FetalTissue_withArtifacts_InteractionStrength.pdf", width = 12, height = 10)
netVisual_circle(cellchat.FT@net$weight, vertex.weight = groupSize.FT, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat.FT@net$weight, vertex.weight = groupSize.FT, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.01)
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks
pdf("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_FetallTissue-individualCelltype.pdf", width = 12, height = 10)
mat.FT <- cellchat.FT@net$weight
# par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat.FT)) {
  mat2 <- matrix(0, nrow = nrow(mat.FT), ncol = ncol(mat.FT), dimnames = dimnames(mat.FT))
  mat2[i, ] <- mat.FT[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.FT, weight.scale = T, edge.weight.max = max(mat.FT), title.name = rownames(mat.FT)[i])
}
dev.off()

# Visualize each signaling pathway
cellchat.FT@netP$pathways

list_available_pathways <- cellchat.FT@netP$pathways

for (pathway in list_available_pathways) {
  print(paste0("=> CCC for pathway: ", pathway))
  
  pdf(paste0("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_FetalTissue_signalingPathway_", pathway, ".pdf"), width = 12, height = 10)
  
  pathways.show <- c(pathway) 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat.FT, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.FT, signaling = pathways.show, layout = "circle")
  
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.FT, signaling = pathways.show, layout = "circle", vertex.label.cex = 0.01)
  
  
  # Chord diagram
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat.FT, signaling = pathways.show, layout = "chord")
  
  # Heatmap
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat.FT, signaling = pathways.show, color.heatmap = "Reds"))
  
  dev.off()
  
}

saveRDS(cellchat.FT, "E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchat/CCC_FetalTissue.RDS")



save.image("E:/Tays_workspace/Demultiplexing_placenta_results/TBannotation_wholeblood/results/cellchate_placenta_TB_annotation_wholeblood.RData")

