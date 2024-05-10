#' scATAC-seq Hematopoiesis data
#'
#' @description
#' single-cell ATAC-seq data from human hematopoiesis (Granja et al. 2019), subset for bone-marrow-derived cells
#'
#' @format  ## `hemato`
#'  A list of three components:
#' (i) counts: a sparse matrix containing the peaks-by-cells counts matrix, unbinarized
#' (ii) meta: cell metadata including cluster and cell-type annotation corresponding to colnames of counts
#' (iii) ranges: a GRanges object corresponding to rownames of counts
#' For testing/example purposes the data were subset to 10000 randomly selected rows and subset for progenitor cells, monocytes, B, T and NK cells.
#'
#' @source <https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds> "hemato"
#'
#' @rdname hemato
"hemato"
