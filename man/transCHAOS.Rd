% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transCHAOS.R
\name{transCHAOS}
\alias{transCHAOS}
\title{Compute transcriptional heterogeneity scores given a scRNA-seq counts matrix and metadata}
\usage{
transCHAOS(
  counts,
  meta,
  colname = colnames(meta)[1],
  n = 100,
  index = NULL,
  subsample = 1
)
}
\arguments{
\item{counts}{A counts matrix representing e.g. a normalised & scaled scRNA-seq matrix}

\item{meta}{A dataframe containing metadata including at least the column on which the data should be grouped, e.g. cluster or cell type.}

\item{colname}{The name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".}

\item{n}{The number of cells to subset for each group/cluster. This defaults to 100.}

\item{index}{The rows/genes on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a vector of numerical indices or a vector or names corresponding to the rownames of interest in "counts".}

\item{subsample}{The number of times cells should be subsampled as "pseudoreplicates" for between-group comparisons. This defaults to 1. If > 1, "n" cells are subsampled from each group the specified number of times.}
}
\value{
A dataframe containing transcriptional heterogeneity scores for each group.
}
\description{
Compute transcriptional heterogeneity scores given a scRNA-seq counts matrix and metadata
}
\examples{
gex <- hemato_rna$data.scale
meta <- hemato_rna$meta
heterogeneity <- transCHAOS(counts=gex, meta=meta, colname="Clusters")
}
