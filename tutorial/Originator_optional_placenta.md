(OPTIONAL) Perform maternal and fetal origin identification in placenta tissues
-------------------------------------------------------------------------------

Since placenta tissues not only contain the mixture of blood and tissue-resident immune cells but also cells from maternal and fetal origins. Orignator provides an additional steps to perform genetic origin identification by incorporating **Freemuxlet** introduced by Popsclde (<https://github.com/statgen/popscle>)

This section can be performed by submitting the job to the Slurm server. This section can be modified based on your local computing environment.

``` bash
#!/bin/bash

#SBATCH --job-name=popscle_PNAS             ## Name of the job for the scheduler
                                       ## make the directive = #SBATCH, not ##SBATCH
#SBATCH --nodes=1                      ## number of nodes you are requesting
#SBATCH --ntasks-per-node=1            ## how many cores do you want to reserve
#SBATCH --time=4-23:30:00              ## Maximum length of time you are reserving the
                                       ## resources for
                                       ## (if job ends sooner, bill is based on time used)
#SBATCH --mem-per-cpu=70g               ## Memory requested per core
#SBATCH --mail-user=qhhuang@umich.edu  ## send email notifications to umich email listed
#SBATCH --mail-type=END                ## when to send email (standard values are:
                                       ## NONE, BEGIN, END, FAIL, REQUEUE, ALL.
                                       ## (See documentation for others)
#SBATCH --output=./%x-%j               ## send output and error info to the file listed
                                       ##(optional: different name format than default)

# I recommend using the following lines to write output to indicate your script is working
if [[ $SLURM_JOB_NODELIST ]] ; then
   echo "Running on"
   scontrol show hostnames $SLURM_JOB_NODELIST
fi

# With SLURM, you can load your modules in the SBATCH script

#  Put your job commands after this line

module load singularity

cd /home/qhhuang/scRNA_PE/Analysis_11_22_2021/

### for i in {1..22} X Y;do echo "${i} chr${i}";done > rename_chrm.txt
### /home/qhhuang/scRNA_PE/Analysis_11_22_2021/bcftools/bcftools annotate filterbyMAF.vcf.gz --r  ename-chrs rename_chrm.txt -Oz -o filterbyMAF_renamed.vcf.gz

###PNAS PE1

for ID in PE2 PE3 PE4
do
  ### sort vcf and bam file order
  /home/qhhuang/scRNA_PE/Analysis_11_22_2021/popscle_helper_tools/sort_vcf_same_as_bam.sh \
  /home/yhdu/hongkong/ALIGNMENT/${ID}/${ID}.bam \
  /home/qhhuang/scRNA_PE/Analysis_11_22_2021/filterbyMAF_renamed.vcf.gz \
  v \
    > /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}.sorted_as_in_bam.vcf

  cp /home/yhdu/hongkong/ALIGNMENT/${ID}/filtered_gene_bc_matrices/hg19/barcodes.tsv /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}_barcodes.tsv
  ### Filter bam for pile up
  /home/qhhuang/scRNA_PE/Analysis_11_22_2021/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
      /home/yhdu/hongkong/ALIGNMENT/${ID}/${ID}.bam \
      /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}_barcodes.tsv \
      /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}.sorted_as_in_bam.vcf \
      /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}_to_demultiplex.filter_bam_file_for_popscle_dsc_pileup.bam

  # Use filtered BAM file for dsc-pileup.
  singularity run popscle_latest.sif "dsc-pileup \
      --sam /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}_to_demultiplex.filter_bam_file_for_popscle_dsc_pileup.bam \
      --vcf /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}.sorted_as_in_bam.vcf \
      --group-list /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}_barcodes.tsv \
      --out /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/${ID}_to_demultiplex.pileup"

  # freemuxlet
  singularity run popscle_latest.sif "freemuxlet \
      --plp /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/${ID}_to_demultiplex.pileup --nsample 2 --out /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/${ID}_freemuxlet.pooled"

  mv /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/${ID}_freemuxlet.pooled.clust1.samples.gz /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/${ID}_freemuxlet.pooled.clust1.samples.gz
  rm /home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/${ID}*.*
  rm /home/qhhuang/scRNA_PE/Analysis_11_22_2021/${ID}*.*
done[thatchau@garmire-gpu01 demultiplexing_project]
```

### Perform data integration and include identified genetic origins into the data.

``` r
##location for downloaded data /home/yhdu/hongkong/ALIGNMENT/
#####load data from own ######
##remove genes that present in less than 3 cells and remove cell with less than 200 features expressed
Ctrl1.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Control1/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_Ctrl1 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/control1_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_Ctrl1 <- as.data.frame(freemuxlet_Ctrl1)
freemuxlet_Ctrl1 <- separate(freemuxlet_Ctrl1, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_Ctrl1) <- freemuxlet_Ctrl1$BARCODE
freemuxlet_Ctrl1 <- freemuxlet_Ctrl1[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_Ctrl1$SNG.BEST.GUESS[freemuxlet_Ctrl1$SNG.BEST.GUESS == 0] <- "Control1_0"
freemuxlet_Ctrl1$SNG.BEST.GUESS[freemuxlet_Ctrl1$SNG.BEST.GUESS == 1] <- "Control1_1"


Ctrl1 <- CreateSeuratObject(counts = Ctrl1.data, min.cells = 3, min.features = 200, project = "Control1")
Ctrl1 <- AddMetaData(Ctrl1, freemuxlet_Ctrl1, col.name = NULL)
Ctrl1$cond <- "Control 1"
Ctrl1$group <- "Controls"
Ctrl1$location <- "Maternal Surface"
Ctrl1$data_source <- "Own"
Ctrl1 <- subset(Ctrl1, subset = DROPLET.TYPE != "DBL")


Ctrl2.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Control2/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_Ctrl2 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/control2_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_Ctrl2 <- as.data.frame(freemuxlet_Ctrl2)
freemuxlet_Ctrl2 <- separate(freemuxlet_Ctrl2, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_Ctrl2) <- freemuxlet_Ctrl2$BARCODE
freemuxlet_Ctrl2 <- freemuxlet_Ctrl2[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_Ctrl2$SNG.BEST.GUESS[freemuxlet_Ctrl2$SNG.BEST.GUESS == 0] <- "Control2_0"
freemuxlet_Ctrl2$SNG.BEST.GUESS[freemuxlet_Ctrl2$SNG.BEST.GUESS == 1] <- "Control2_1"

Ctrl2 <- CreateSeuratObject(counts = Ctrl2.data, min.cells = 3, min.features = 200, project = "Control2")
Ctrl2 <- AddMetaData(Ctrl2, freemuxlet_Ctrl2, col.name = NULL)
Ctrl2$cond <- "Control 2"
Ctrl2$group <- "Controls"
Ctrl2$location <- "Fetal Surface"
Ctrl2$data_source <- "Own"
Ctrl2 <- subset(Ctrl2, subset = DROPLET.TYPE != "DBL")

Ctrl3.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Control3/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_Ctrl3 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/control3_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_Ctrl3 <- as.data.frame(freemuxlet_Ctrl3)
freemuxlet_Ctrl3 <- separate(freemuxlet_Ctrl3, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_Ctrl3) <- freemuxlet_Ctrl3$BARCODE
freemuxlet_Ctrl3 <- freemuxlet_Ctrl3[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_Ctrl3$SNG.BEST.GUESS[freemuxlet_Ctrl3$SNG.BEST.GUESS == 0] <- "Control3_0"
freemuxlet_Ctrl3$SNG.BEST.GUESS[freemuxlet_Ctrl3$SNG.BEST.GUESS == 1] <- "Control3_1"

Ctrl3 <- CreateSeuratObject(counts = Ctrl3.data, min.cells = 3, min.features = 200, project = "Control3")
Ctrl3 <- AddMetaData(Ctrl3, freemuxlet_Ctrl3, col.name = NULL)
Ctrl3$cond <- "Control 3"
Ctrl3$group <- "Controls"
Ctrl3$location <- "Intermediate Surface"
Ctrl3$data_source <- "Own"
Ctrl3 <- subset(Ctrl3, subset = DROPLET.TYPE != "DBL")

Case1.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Garmire1/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_case1 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/case1_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_case1 <- as.data.frame(freemuxlet_case1)
freemuxlet_case1 <- separate(freemuxlet_case1, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_case1) <- freemuxlet_case1$BARCODE
freemuxlet_case1 <- freemuxlet_case1[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_case1$SNG.BEST.GUESS[freemuxlet_case1$SNG.BEST.GUESS == 0] <- "case1_0"
freemuxlet_case1$SNG.BEST.GUESS[freemuxlet_case1$SNG.BEST.GUESS == 1] <- "case1_1"

Case1 <- CreateSeuratObject(counts = Case1.data, min.cells = 3, min.features = 200, project = "Case1")
Case1 <- AddMetaData(Case1, freemuxlet_case1, col.name = NULL)
Case1$cond <- "Case 1"
Case1$group <- "Cases"
Case1$location <- "Maternal Surface"
Case1$data_source <- "Own"
Case1 <- subset(Case1, subset = DROPLET.TYPE != "DBL")


Case2.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Garmire2/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_case2 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/case2_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_case2 <- as.data.frame(freemuxlet_case2)
freemuxlet_case2 <- separate(freemuxlet_case2, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_case2) <- freemuxlet_case2$BARCODE
freemuxlet_case2 <- freemuxlet_case2[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_case2$SNG.BEST.GUESS[freemuxlet_case2$SNG.BEST.GUESS == 0] <- "case2_0"
freemuxlet_case2$SNG.BEST.GUESS[freemuxlet_case2$SNG.BEST.GUESS == 1] <- "case2_1"

Case2 <- CreateSeuratObject(counts = Case2.data, min.cells = 3, min.features = 200, project = "Case2")
Case2 <- AddMetaData(Case2, freemuxlet_case2, col.name = NULL)
Case2$cond <- "Case 2"
Case2$group <- "Cases"
Case2$location <- "Fetal Surface"
Case2$data_source <- "Own"
Case2 <- subset(Case2, subset = DROPLET.TYPE != "DBL")

Case3.data <- Read10X(data.dir = "/home/qhhuang/scRNA_PE/scRNA_data/Sample_Garmire3/outs/filtered_gene_bc_matrices/hg19")
#load freemuxlet result to remove doublet
freemuxlet_case3 <- fread(paste0("zcat < ", "/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/case3_freemuxlet.pooled.clust1.samples.gz"));
freemuxlet_case3 <- as.data.frame(freemuxlet_case3)
freemuxlet_case3 <- separate(freemuxlet_case3, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
row.names(freemuxlet_case3) <- freemuxlet_case3$BARCODE
freemuxlet_case3 <- freemuxlet_case3[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
freemuxlet_case3$SNG.BEST.GUESS[freemuxlet_case3$SNG.BEST.GUESS == 0] <- "case3_0"
freemuxlet_case3$SNG.BEST.GUESS[freemuxlet_case3$SNG.BEST.GUESS == 1] <- "case3_1"

Case3 <- CreateSeuratObject(counts = Case3.data, min.cells = 3, min.features = 200, project = "Case3")
Case3 <- AddMetaData(Case3, freemuxlet_case3, col.name = NULL)
Case3$cond <- "Case 3"
Case3$group <- "Cases"
Case3$location <- "Intermediate Surface"
Case3$data_source <- "Own"
Case3 <- subset(Case3, subset = DROPLET.TYPE != "DBL") 

data_list <- list(Ctrl1, Ctrl2, Ctrl3, Case1, Case2, Case3)

#####load data from PNAS ######
data_name <- c("PE1", "PE2", "PE3", "PE4", "PN1", "PN2", "PN3C", "PN4C")

for (dn in data_name){
  file_path1 = paste("/home/yhdu/hongkong/ALIGNMENT/", dn, "/filtered_gene_bc_matrices/hg19", sep = "")
  temp.data1 <- Read10X(data.dir = file_path1)
  file_path2 = paste("/home/qhhuang/scRNA_PE/Analysis_11_22_2021/result/essential_result/", dn, "_freemuxlet.pooled.clust1.samples.gz", sep = "")
  temp.result1 = fread(paste0("zcat < ", file_path2))
  temp.result1 <- as.data.frame(temp.result1)
  temp.result1 <- separate(temp.result1, BARCODE, sep = "-", c('BARCODE', 'BARCODE-2'))
  row.names(temp.result1) <- temp.result1$BARCODE
  temp.result1 <- temp.result1[, c("BARCODE", "SNG.BEST.GUESS", "DROPLET.TYPE")]
  temp.result1$SNG.BEST.GUESS[temp.result1$SNG.BEST.GUESS == 0] <- paste0(dn, "_0")
  temp.result1$SNG.BEST.GUESS[temp.result1$SNG.BEST.GUESS == 1] <- paste0(dn, "_1")
  
  temp.data2 <- CreateSeuratObject(counts = temp.data1, min.cells = 3, min.features = 200, project = dn)
  temp.data2 <- AddMetaData(temp.data2, temp.result1, col.name = NULL)
  temp.data2$location <- "Intermediate Surface"
  if (dn %in% c("PE1", "PE2", "PE3", "PE4")){
    temp.data2$cond <- dn
    temp.data2$group <- "Cases"
  }else{
    temp.data2$cond <- dn
    temp.data2$group <- "Controls"
  }
  temp.data2$data_source <- "PNAS"
  temp.data2 <- subset(temp.data2, subset = DROPLET.TYPE != "DBL")
  #append datalist
  data_list <- c(data_list, temp.data2)
}

#QC remove outliers based on QC metrics
for(j in c(1:length(data_list))){
  pdf(paste('QC_plots_', j, '.pdf', sep = ''), width = 12, height = 5)
  data_list[[j]][["percent.mt"]] <- PercentageFeatureSet(data_list[[j]], pattern = "^MT-")
  print(VlnPlot(data_list[[j]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  plot1 <- FeatureScatter(data_list[[j]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data_list[[j]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  dev.off()
}

##NOTE PNAS data do not have mito genes

#Filter mito < 10% and Count > 200 and < 6000 based on the visual inspection of graph
#Standard pre-processing
for(j in c(1:length(data_list))){
  data_list[[j]] <- subset(data_list[[j]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  if(j %in% c(1:6)){
    data_list[[j]] <- SCTransform(data_list[[j]], vars.to.regress = "percent.mt", verbose = FALSE)
  }else{
    data_list[[j]] <- SCTransform(data_list[[j]], verbose = FALSE)
  }
}

save(data_list, file = "/home/qhhuang/scRNA_PE/Analysis_4_25_2022/data_list.RData")

###integrate###
PL.features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
assay_select <- rep("SCT", 14)
PL.anchors <- FindIntegrationAnchors(object.list = data_list, dims = 1:30, anchor.features = PL.features, assay = assay_select)
PL.integrated <- IntegrateData(anchorset = PL.anchors, dims = 1:30)

###Further analysis on the integrated assays
DefaultAssay(PL.integrated) <- "integrated"
PL.integrated <- ScaleData(PL.integrated, verbose = FALSE)
PL.integrated <- RunPCA(PL.integrated, verbose = F)
# UMAP and Clustering
PL.integrated <- RunUMAP(PL.integrated, reduction = "pca", dims = 1:30)
PL.integrated <- FindNeighbors(PL.integrated, reduction = "pca", dims = 1:30)
# PL.integrated <- FindClusters(PL.integrated, resolution = 0.4)
PL.integrated <- FindClusters(PL.integrated, resolution = 0.6)

pdf('Visualization of the clustering3.pdf', width = 12, height = 6)
DimPlot(PL.integrated, reduction = "umap", group.by = "integrated_snn_res.0.6", split.by = "data_source", label = TRUE)
dev.off()

# pdf('Visualization of the clustering.pdf', width = 10, height = 6)
# DimPlot(PL.integrated, reduction = "umap", split.by = "data_source", label = TRUE)
# dev.off()

###annotate cell-type
DefaultAssay(PL.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
PL.integrated <- NormalizeData(PL.integrated, verbose = FALSE)

pdf("Gene Marker for VCT.pdf")
FeaturePlot(object = PL.integrated, features = c("PAGE4", "PEG10", "ISYNA1"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for SCT.pdf")
FeaturePlot(object = PL.integrated, features = c("CGA", "CYP19A1"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for EVT.pdf")
FeaturePlot(object = PL.integrated, features = c("HLA-G",   "HTRA4"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for FB1.pdf")
FeaturePlot(object = PL.integrated, features = c("DLK1",    "EGFL6", "GPC3", "GPC3"),min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for FB2.pdf")
FeaturePlot(object = PL.integrated, features = c("DKK1","IGFBP5", "IGFBP1"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

# cluster_22_markers <- FindMarkers(PL.integrated, ident.1 = 22 , ident.2 = NULL, only.pos = TRUE)
# head(cluster_22_markers)

pdf("Gene Marker for Endothelial.pdf")
FeaturePlot(object = PL.integrated, features = c("CD34", "VWF", "CLEC14A",  "ECSCR"),min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for EB.pdf")
FeaturePlot(object = PL.integrated, features = c("HBB", "HBG1", "HBA1"),min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for Hofbauer cells.pdf")
FeaturePlot(object = PL.integrated, features = c("LYVE1", "DAB2", "CCL4", "CSF1R", "CD163", "CD209"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for DCs.pdf")
FeaturePlot(object = PL.integrated, features = c("FCER1A", "CST3"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for T cells.pdf")
FeaturePlot(object = PL.integrated, features = c("IL7R", "CCR7"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for NK cells.pdf")
FeaturePlot(object = PL.integrated, features = c("GNLY", "NKG7"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for FCGR3A+ Mono cells.pdf")
FeaturePlot(object = PL.integrated, features = c("FCGR3A", "MS4A7"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for CD14+ Mono cells.pdf")
FeaturePlot(object = PL.integrated, features = c("CD14", "LYZ"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

pdf("Gene Marker for B.pdf")
FeaturePlot(object = PL.integrated, features = c("MS4A1"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.5, combine = F)
dev.off()

###change annotation villous trophoblast to cytotrophoblasts

new.ident <- c("Monocyte",  "Extravillous Trophoblast", "Fibroblast Type 1", "Cytotrophoblast", "Fibroblast Type 2", "Fibroblast Type 1", "Extravillous Trophoblast", "Macrophage (HB)", "Cytotrophoblast", "Cytotrophoblast", "Monocyte", "Syncytiotrophoblast", "Vascular Endothelial", "Erythrocyte", "T-cells", "Fibroblast Type 1", "Extravillous Trophoblast", "Fibroblast Type 2", "Cytotrophoblast", "Fibroblast Type 1", "Syncytiotrophoblast", "Fibroblast Type 1", "Macrophage (HB)", "Natural Killer cells", "Fibroblast Type 1", "noise")

names(new.ident) <- levels(x = PL.integrated)
PL.integrated <- RenameIdents(object = PL.integrated, new.ident)
PL.integrated$annotation <- Idents(PL.integrated)
save(PL.integrated, file = "/home/qhhuang/scRNA_PE/Analysis_4_25_2022/PL.integrated_4_25_2022.RData")

###cluster to remove from downstream analysis: 25

PL.integrated_1 <- subset(PL.integrated, subset = annotation != "noise")
PL.integrated_2 <- subset(PL.integrated_1, subset = DROPLET.TYPE == "SNG")

save(PL.integrated_2, file = "/home/qhhuang/scRNA_PE/Analysis_4_25_2022/PL.integrated_SNG_4_25_2022.RData")
```

### Maternal and fetal origin

Based on the cell type information and the results from Freemuxlet, maternal and fetal origin identification can be performed as follows.

``` r
##Assign fetal and maternal origin based on trophoblast cluster assignment
pdf('Visualization of the SNP prediction.pdf', width = 30, height = 25)
DimPlot(PL.integrated_2, reduction = "umap",  group.by = "annotation", split.by = "SNG.BEST.GUESS", label = TRUE)
dev.off()

PL_meta <- PL.integrated_2[[]]
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case1_1"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case1_0"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case2_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case2_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case3_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "case3_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Control1_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Control2_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Control2_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Control3_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Control3_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE1_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE1_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE2_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE2_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE3_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE3_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE4_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PE4_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN1_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN1_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN2_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN2_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN3C_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN3C_1"] <- "Maternal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN4C_0"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "PN4C_1"] <- "Maternal"
##Assign some scattered trophoblasts in maternal origin to fetal origin exclusively
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Maternal" & PL_meta$annotation == "Extravillous Trophoblast"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Maternal" & PL_meta$annotation == "Cytotrophoblast"] <- "Fetal"
PL_meta$SNG.BEST.GUESS[PL_meta$SNG.BEST.GUESS == "Maternal" & PL_meta$annotation == "Syncytiotrophoblast"] <- "Fetal"

PL.integrated_2[["MF_Origin"]] <- PL_meta$SNG.BEST.GUESS

###compartment assignment
PL_meta$annotation <- as.character(PL_meta$annotation)
PL_meta$compartment <- PL_meta$annotation
PL_meta$compartment[PL_meta$annotation == "Monocyte" | PL_meta$annotation == "Macrophage (HB)" | PL_meta$annotation == "Erythrocyte" | PL_meta$annotation == "T-cells" | PL_meta$annotation == "Natural Killer cells"] <- "Immune"
PL_meta$compartment[PL_meta$annotation == "Vascular Endothelial" | PL_meta$annotation == "Fibroblast Type 1" | PL_meta$annotation == "Fibroblast Type 2"] <- "Stromal"
PL_meta$compartment[PL_meta$annotation == "Extravillous Trophoblast" | PL_meta$annotation == "Syncytiotrophoblast" | PL_meta$annotation == "Cytotrophoblast"] <- "Trophoblast"
PL.integrated_2[["compartment"]] <- PL_meta$compartment

save(PL.integrated_2, file = "/home/qhhuang/scRNA_PE/Analysis_4_25_2022/PL.integrated_SNG_4_25_2022.RData")

pdf("/home/qhhuang/scRNA_PE/Analysis_4_25_2022/MF_Origin_annotation.pdf", width = 6, height = 5)
DimPlot(PL.integrated_2, group = "MF_Origin", cols = c('Fetal' = 'dodgerblue4', "Maternal" = 'grey'), reduction = "umap")
dev.off()

pdf("/home/qhhuang/scRNA_PE/Analysis_4_25_2022/Celltypes_annotation.pdf", width = 7, height = 5)
DimPlot(PL.integrated_2, group = "annotation", reduction = "umap", label=TRUE, repel = TRUE)
dev.off()

pdf('Visualization of the Batch Integration.pdf', width = 13, height = 6)
DimPlot(PL.integrated_2, reduction = "umap", group.by = "annotation", split.by = "data_source", label = TRUE, repel = TRUE)
dev.off()
```

The output from here will be further used to perform blood and tissue-immune cell identification.
