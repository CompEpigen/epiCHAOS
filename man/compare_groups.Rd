% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential.R
\name{compare_groups}
\alias{compare_groups}
\title{Function to compare epiCHAOS scores between two groups with permutation test to compute p-value}
\usage{
compare_groups(group1, group2, niter = 1000)
}
\arguments{
\item{group1}{A binarised peaks-by-cells matrix (or other single cell epigenomics matrix) for the first comparison group.}

\item{group2}{A binarised peaks-by-cells matrix (or other single cell epigenomics matrix) for the second comparison group.}

\item{niter}{A numeric value indicating the number of permutations. Defaults to 1000.}
}
\value{
A list containing (i) "result.mix": epiCHAOS scores from the permuted mixtures, and
(ii) "result.test": epiCHAOS scores from the test groups.
(iii) "pval": a p-value for the permutation test.
}
\description{
Function to compare epiCHAOS scores between two groups with permutation test to compute p-value
}
