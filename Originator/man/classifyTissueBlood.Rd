% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Originator.R
\name{classifyTissueBlood}
\alias{classifyTissueBlood}
\title{Perform blood and tissue-resident immune cell identification using Originator}
\usage{
classifyTissueBlood(
  data,
  celltype,
  output_dir,
  use_ident = "orig.ident",
  use_celltype = "unified_celltype",
  into_col = "origin_tb",
  plot = T,
  offset = -0.1,
  reduction = "umap"
)
}
\arguments{
\item{data}{Data}

\item{output_dir}{Output path to save all results}

\item{use_ident}{Column name indicating query or reference}

\item{use_celltype}{Column name indicating cell type}

\item{into_col}{Column to save classification result}

\item{plot}{Generates intermediate plots and save in output path}

\item{offset}{Offset of embedding borders}

\item{celtype}{Cell type of interest}
}
\value{
Tissue blood origin annotated data
}
\description{
Perform blood and tissue-resident immune cell identification using Originator
}
