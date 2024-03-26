# Originator: Computational Framework Separating Single-Cell RNA-Seq by Genetic and Contextual Information

Single-cell RNA sequencing (scRNA-Seq) data from tissues are prone to blood contamination in sample preparation. Moreover, some tissue samples comprise cells of different genetic makeups. These issues require rigorous preprocessing and cell filtering prior to the downstream functional analysis. We propose a new computational framework, Originator, which deciphers single cells into different genetic origins and separates blood cells from tissue-resident cells in the scRNA-Seq data. We show that this pipeline improves downstream data analysis, exemplified by the pancreatic ductal adenocarcinoma and placenta tissues.

<img src="./image/originator_pipeline.png" alt="originator_pipeline" width="500"/>

# Installation
Install dependencies
``R
install.packages("tidyr")
install.packages('Seurat')
install.packages("RColorBrewer")
install.packages("matrixStats")
```
Install Originator
``R
remotes::install_github("thatchayut/Originator/Originator/")
```
