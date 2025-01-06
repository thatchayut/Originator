setwd("/nfs/dcmb-lgarmire/johntao/originator")
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
load('PL.integrated_SNG_4_25_2022.RData')
PN3C<-subset(PL.integrated_2,subset=orig.ident=="PN3C")
file_path2 = paste("/nfs/dcmb-lgarmire/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/", "PN3C", "_freemuxlet.pooled.clust1.samples.gz", sep = "")
temp.result1 = fread(paste0("zcat < ", file_path2))
temp.result1 <- as.data.frame(temp.result1)
temp.result1 <- separate(temp.result1, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(temp.result1) <- temp.result1$BARCODE
temp.result1 <- temp.result1[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
temp.result1$SNG.BEST.GUESS[temp.result1$SNG.BEST.GUESS == 0] <- paste0("PN3C", "_0")
temp.result1$SNG.BEST.GUESS[temp.result1$SNG.BEST.GUESS == 1] <- paste0("PN3C", "_1")

temp.result1.subset<-subset(temp.result1,subset=BARCODE %in% PN3C@meta.data$BARCODE)
temp.result1.subset$BARCODE<-paste0(temp.result1.subset$BARCODE,"_13")

temp.result1.subset<-temp.result1.subset[match(row.names(PN3C@meta.data), temp.result1.subset$BARCODE),]

PN3C_meta<-PN3C@meta.data
PN3C_meta<-cbind(PN3C_meta,temp.result1.subset)

PN3C_meta$SNG.BEST.GUESS.MF<-rep("NA",nrow(PN3C_meta))
PN3C_meta$SNG.BEST.GUESS.MF[PN3C_meta$SNG.BEST.GUESS == "PN3C_0"] <- "Fetal"
PN3C_meta$SNG.BEST.GUESS.MF[PN3C_meta$SNG.BEST.GUESS == "PN3C_1"] <- "Maternal"


PN3C@meta.data<-PN3C_meta
saveRDS(PN3C,"PN3C_withSNGBestGuess.RDS")

PN3C_trophoblast<-subset(PN3C,subset=compartment=="Trophoblast")

scSplit<-read.csv("./PN3CscSplit/scSplit_result.csv")
result<-scSplit %>%
  separate(Barcode.Cluster,into=c("Barcode","Group"),sep='\t')
result$Barcode<-gsub("-1","_13",result$Barcode)

PN3C_trophoblast$barcode_dummy <- paste0(PN3C_trophoblast$BARCODE, "_13")

PN3C_trophoblast_subset<-subset(PN3C_trophoblast,subset = barcode_dummy %in% result$Barcode)

result<-result[match(PN3C_trophoblast_subset@meta.data$barcode_dummy,result$Barcode),]
PN3C_trophoblast_subset_meta<-PN3C_trophoblast_subset@meta.data


PN3C_trophoblast_subset_meta<-cbind(PN3C_trophoblast_subset_meta,result)

PN3C_trophoblast_subset@meta.data<-PN3C_trophoblast_subset_meta

#### Assigning class for scSplit
test <- subset(PN3C_trophoblast_subset, subset = Group != "DBL-1")

# Use table(test$Group) to observe number of cells in each scSplit class
##  assign the larger group as fetal side
## Make sure to save R object using saveRDS()

