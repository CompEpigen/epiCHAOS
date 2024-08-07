
#' @importFrom magrittr %>%
magrittr::`%>%`

#' @importFrom IRanges subsetByOverlaps
IRanges::subsetByOverlaps

#' @importFrom GenomicRanges GRanges
GenomicRanges::GRanges


#' Compare epiCHAOS scores across a set of genomic regions
#'
#' @param counts a binarised sparse matrix containing regions/peaks-by-cells scATAC or other single-cell epigenomics dataset. Rownames should contain the genomic positions of peaks or regions in the format: ("chrX:start-end")
#' @param regions a GRangesList object containing the ranges corresponding to regions of interest, e.g. the ENCODE TFBS or any other custom regions dataset
#' @param n the number of peaks to be subset for each region type
#' @param min.row the minimum number of rows/peaks for a region type. If the number of rows overlapping with regions os a specified type is too few, epiCHAOS scores are not computed for that region type.
#'
#' @return a dataframe containing epiCHAOS scores for each selected region in the input data
#' @export
#'
#' @examples
#' counts <- hemato_atac$counts
#' rownames(counts) <- paste0(hemato_atac$ranges)
#' regions <- encode_tfbs$ranges
#' per.region.scores <- compare_regions(counts=counts, regions=regions, n=100, min.row=20)
compare_regions <- function(counts, regions, n=2000, min.row=20) {

  matrices <- list()

  #--- Genomic Ranges for the input data
  data.gr <- rownames(counts) %>% GRanges()
  names(data.gr) <- rownames(counts)

  regions.interest <- names(regions)

  #--- loop through each region
  for (region in regions.interest) {

    #--- get GRanges object containing selected regions and subset the peak/counts matrix for the selected regions
    temp <- regions[[region]]
    select.sites <- IRanges::subsetByOverlaps(data.gr, temp) %>% names()
    temp <- counts[select.sites, ]

    #--- remove region types with small number of peaks
    if (nrow(temp)<min.row) { message("too few overlapping peaks, moving to next region"); next }

    #--- subset to n regions if specified
    if (nrow(temp)<n) { temp <- temp } else { temp <- temp[sample(nrow(temp), n, replace = F), ]}

    #--- binarise if not already
    temp[temp>1] <- 1
    matrices[[region]] <- temp

  }

  #--- compute epiCHAOS scores
  het <- compute_eITH(matrices)
  return(het)

}

