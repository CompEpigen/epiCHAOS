


#' @importFrom magrittr %>%
magrittr::`%>%`

#' @importFrom IRanges subsetByOverlaps
IRanges::subsetByOverlaps

#' @importFrom GenomicRanges GRanges
GenomicRanges::GRanges

#' @import TxDb.Hsapiens.UCSC.hg19.knownGene

#' @import TxDb.Hsapiens.UCSC.hg38.knownGene

#' @importFrom msigdbr msigdbr
msigdbr::msigdbr

#' @importFrom GenomicFeatures promoters genes
GenomicFeatures::promoters
GenomicFeatures::genes


#' Compare epiCHAOS scores across selected gene sets or biological pathways
#'
#' @param counts a binarised sparse matrix containing regions/peaks-by-cells scATAC or other single-cell epigenomics dataset. Rownames should contain the genomic positions of peaks or regions in the format: "chrX:start-end"
#' @param min.row the minimum number of rows/peaks overlapping with promoters of selected gene sets. If the number of rows overlapping with promoters of a selected gene set is too few, epiCHAOS scores are not computed for that gene set.
#' @param set.names a vector of gene sets from msigdbr on which to compare epiCHAOS scores
#' @param genome the genome build e.g. "hg19", "hg38" on which to define promoter-overlapping regions. This defaults to hg19
#' @param tss.span a numeric vector containing the number of bp upstream and downstream of the TSS by which promoters should be defined. This defaults to c(1500, 500) i.e. 1500bp upstream and 500bp downstream of the TSS.
#'
#' @return a dataframe containing epiCHAOS scores for each selected gene set
#' @export
#'
#' @examples
#' library(magrittr)
#' counts <- hemato_atac$counts
#' rownames(counts) <- paste0(hemato_atac$ranges)
#'
#' # load hallmark gene ontologies for human from msigdb
#' gene_sets = msigdbr::msigdbr(species = "Homo sapiens")
#' gene_sets <- gene_sets[gene_sets$gs_cat=="H",]
#' set.names <- gene_sets$gs_name %>% unique()
#' het <- compare_pathways(counts=counts, set.names=set.names, min.row=5,
#' genome="hg19", tss.span=c(1500, 500))
compare_pathways <- function(counts, set.names, min.row=10, genome="hg19", tss.span=c(1500, 500)) {

  #--- get gene sets from msigdbr, subset for gene sets of interest
  gene_sets <- msigdbr::msigdbr(species = "Homo sapiens")
  gene_sets <- gene_sets[gene_sets$gs_name %in% set.names,]


  #--- select hg19/hg38 promoters
  if (genome=="hg19") { txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene } else if (genome=="hg38") { txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene }

  promoters <- GenomicFeatures::promoters(genes(txdb), upstream = tss.span[1], downstream = tss.span[2])

  #--- to hold peaks matrices for each selected gene set/pathway
  matrices <- list()

  #--- binarise counts if not already & create GRanges object from rownames
  counts[counts>1] <- 1
  ranges <- rownames(counts) %>% GRanges()
  names(ranges) <- rownames(counts)

  #--- loop through gene sets
  message("computing epiCHAOS scores for selected gene sets/pathways")
  for (set in set.names) {

    #print(set)

    #--- select promoters for gene set of interest
    genes <- gene_sets$entrez_gene[gene_sets$gs_name==set] %>% intersect(names(promoters))
    select.promoters <- promoters[genes,]
    select.promoters <- IRanges::subsetByOverlaps(ranges, select.promoters) %>% names()

    #--- exclude gene sets with small number of associated genes/peaks
    if (length(select.promoters)<min.row) { next }

    matrices[[set]] <- counts[select.promoters,]

  }

  #--- compute epiCHAOS scores
  het <- compute_eITH(matrices)

  return(het)
}



