% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_regions.R
\name{compare_regions}
\alias{compare_regions}
\title{Compare epiCHAOS scores across a set of genomic regions}
\usage{
compare_regions(counts, regions, n = 2000, min.row = 20)
}
\arguments{
\item{counts}{a binarised sparse matrix containing regions/peaks-by-cells scATAC or other single-cell epigenomics dataset. Rownames should contain the genomic positions of peaks or regions in the format: ("chrX:start-end")}

\item{regions}{a GRangesList object containing the ranges corresponding to regions of interest, e.g. the ENCODE TFBS or any other custom regions dataset}

\item{n}{the number of peaks to be subset for each region type}

\item{min.row}{the minimum number of rows/peaks for a region type. If the number of rows overlapping with regions os a specified type is too few, epiCHAOS scores are not computed for that region type.}
}
\value{
a dataframe containing epiCHAOS scores for each selected region in the input data
}
\description{
Compare epiCHAOS scores across a set of genomic regions
}
\examples{
counts <- hemato_atac$counts
rownames(counts) <- paste0(hemato_atac$ranges)
regions <- encode_tfbs$ranges
per.region.scores <- compare_regions(counts=counts, regions=regions, n=100, min.row=20)
}
