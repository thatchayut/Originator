# Originator: Computational Framework Separating Single-Cell RNA-Seq by Genetic and Contextual Information

Single-cell RNA sequencing (scRNA-Seq) data from tissues are prone to blood contamination in sample preparation. Moreover, some tissue samples comprise cells of different genetic makeups. These issues require rigorous preprocessing and cell filtering prior to the downstream functional analysis. We propose a new computational framework, Originator, which deciphers single cells into different genetic origins and separates blood cells from tissue-resident cells in the scRNA-Seq data. 

![alt text](image/originator_pipeline.png)
**Illustration of Originator framework** The input data are the scRNA-Seq experiment on tissue sections. (1) data preprocessing and cell type annotation. (2) separating barcoded cells into different origins by blood vs tissue residents context and optionally by inferred genotype information. (3) using the results in steps (1) and (2) to dissect tissue heterogeneity. (4) the functional downstream analyses with respect to cells’ origins.

Please check **[our latest preprint](https://www.biorxiv.org/content/10.1101/2024.04.04.588144v3)**

# Installation

## Install using `renv`

We recommend installing Originator using [renv](https://rstudio.github.io/renv/articles/renv.html), an R package virtual environment manager that simplifies installing packages from various sources and maintaining consistent dependencies.

Install `renv` package manager and activate the environment

```r
install.packages("renv")
renv::init()
renv::activate()
```

Restart your R session and install Originator

```r
renv::install("lanagarmire/Originator/Originator/")
```

## Install using `devtools`

Alternatively, install Originator using `devtools` or `remotes`. Make sure Seurat V4 is installed before proceeding, you may follow the instructions [here](https://satijalab.org/seurat/articles/install_v5.html#install-seurat-v4).

```r
install.packages("devtools")
devtools::install_github("lanagarmire/Originator/Originator/")
# Or use remotes
# install.packages("remotes")
# remotes::install_github("lanagarmire/Originator/Originator/")
```

# Tutorial
- [Perform blood and tissue-resident immune cell identification using Originator](tutorial/Originator_tutorial.md)
- [(OPTIONAL) Perform maternal and fetal origin identification in placenta tissues](tutorial/Originator_optional_placenta.md)


## Contributors
- [Thatchayut Unjitwattana](https://github.com/thatchayut) (thatchau@umich.edu)
- Qianhui Huang (qhhuang@umich.edu)
- [Yiwen Yang](https://github.com/remisiki) (yangiwen@umich.edu)
- Leyang Tao (johntao@umich.edu)
- [Youqi Yang](https://github.com/youqiy) (youqi@umich.edu)
- [Mengtian Zhou](https://github.com/mengtzhou) (mtzhouzz@umich.edu)
- [Yuheng Du](https://github.com/yhdu36) (yhdu@umich.edu)
- [Lana Garmire](https://github.com/lanagarmire) (http://garmiregroup.org/)
