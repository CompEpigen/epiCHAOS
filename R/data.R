#' Example scATAC-seq data
#'
#' @description
#' single-cell ATAC-seq data from human hematopoiesis (Granja et al. 2019), subset for bone-marrow-derived cells, subset for 5000 peaks.
#'
#' @format  ## `hemato_atac`
#'  A list of three components:
#' (i) counts: a sparse matrix containing the peaks-by-cells counts matrix, unbinarized, subset for 5000 peaks
#' (ii) meta: cell metadata including cluster and cell-type annotation corresponding to colnames of counts
#' (iii) ranges: a GRanges object corresponding to rownames of counts
#' For testing/example purposes the data were subset to 10000 randomly selected rows and subset for 100 progenitor cells, B, T and NK cells.
#'
#' @source <https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds> "hemato_atac"
#'
#' @rdname hemato_atac
"hemato_atac"

#' Example scRNA-seq data
#'
#' @description
#' single-cell RNA-seq data from human hematopoiesis (Granja et al. 2019), subset for bone-marrow-derived cells, subset for 5000 genes to use for example purposes
#'
#' @format  ## `hemato_rna`
#'  A list of two components:
#' (i) data.scale: a sparse matrix containing scaled normalised scRNA-seq counts for subset of 100 HSCs, B-cells and T-cells
#' (ii) meta: cell metadata including cluster and cell-type annotation corresponding to colnames of counts
#'
#' @source <https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scRNA-Healthy-Hematopoiesis-191120.rds> "hemato_rna"
#'
#' @rdname hemato_rna
"hemato_rna"

#' Example regions set
#'
#' @description
#' Genomic coordinates in GRanges format for a range of chromatin/transcription factor binding sites from the ENCODE TFBS project, subset for 10 regions to use as example
#'
#' @format  ## `encode_tfbs`
#'  A list of two components:
#' (i) ranges: a GRangesList object containing a GRanges object corresponding to each narrowPeak file from ENCODE TFBS
#' (ii) anno: a dataframe with four columns; filename (corresponding to names of (i)), celltype, treatment and antibody
#'
#' @source <http://cloud.databio.org/regiondb/> "encode_tfbs"
#'
#' @rdname encode_tfbs
"encode_tfbs"
