
#' @importFrom stats dist na.omit sd
stats::dist
stats::na.omit
stats::sd

#' Calculate transcriptional heterogeneity in single-cell RNA-seq data
#'
#' @param x a list of normalised scRNA-seq matrices, each matrix corresponding to a cell group/cluster of interest
#'
#' @return a dataframe containing transcriptional heterogeneity scores
#' @export
#'
#' @examples
#' gex <- hemato_rna$data.scale
#' meta <- hemato_rna$meta
#' matrices <- create.group.matrices(gex, meta, colname="Clusters", binarise=FALSE)
#' heterogeneity <- compute.tITH(matrices)
compute.tITH <- function(x) {

  #--- create list to hold heterogeneity scores per celltype/condition
  dists <- list()
  for (i in names(x)) {
    message(paste0("computing pairwise distances for ", i))
    temp <-  x[[i]] %>% t() %>% dist() %>% as.matrix() %>% colMeans()

    #--- check for and remove outliers in dists[[i]]
    bound <- mean(temp) + 3*sd(temp)
    dists[[i]] <- temp[temp<bound]

  }

  #--- create dataframe to hold heterogeneity scores
  message("Compiling scores to dataframe")
  het <- data.frame(het=unlist(dists))  # this assumes that all datasets have an equal number of cells
  state <- c()
  for (cond in names(dists)) {
    state <- c(state, rep(cond, length(dists[[cond]])))
  }

  het$state <- state

  for (state in unique(het$state)) {
    het$het[het$state==state] <- mean(het$het[het$state==state])
  }

  het <- het[,c("state", "het")] %>% unique()
  het <- na.omit(het)

  #--- transform to a range of 0-1
  het$het <- (het$het - min(het$het)) / (max(het$het) - min(het$het))

  #--- return heterogeneity score
  return(het)
}


#' Compute transcriptional heterogeneity scores given a scRNA-seq counts matrix and metadata
#'
#' @param counts A counts matrix representing e.g. a normalised & scaled scRNA-seq matrix
#' @param meta A dataframe containing metadata including at least the column on which the data should be grouped, e.g. cluster or cell type.
#' @param colname The name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".
#' @param n The number of cells to subset for each group/cluster. This defaults to 100.
#' @param index The rows/genes on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a vector of numerical indices or a vector or names corresponding to the rownames of interest in "counts".
#' @param subsample The number of times cells should be subsampled as "pseudoreplicates" for between-group comparisons. This defaults to 1. If > 1, "n" cells are subsampled from each group the specified number of times.
#'
#' @return A dataframe containing transcriptional heterogeneity scores for each group.
#' @export
#'
#' @examples
#' gex <- hemato_rna$data.scale
#' meta <- hemato_rna$meta
#' heterogeneity <- transCHAOS(counts=gex, meta=meta, colname="Clusters")
transCHAOS <- function(counts, meta, colname=colnames(meta)[1], n=100, index=NULL, subsample=1) {

  # bind variables
  state <- het <- NULL

  #--- create per-group matrices
  message("creating group matrices")

  #--- if subsample is kept at 1, epiCHAOS scores are computed once per group
  if (subsample==1) {
    matrices <- create.group.matrices(counts = counts, meta=meta, colname = colname, n = n, index = index)

    #--- otherwise, a specified number of subsamples of cells are taken for calculation as "pseudoreplicates"
  } else {

    matrices <- list()
    for (rep in 1:subsample) {
      set.seed(rep)
      rep.matrices <- create.group.matrices(counts = counts, meta = meta, colname = colname, n = n, index = index)
      names(rep.matrices) <- paste0(names(rep.matrices), "-", rep)
      matrices <- c(matrices, rep.matrices)
    }
  }

  #--- compute epiCHAOS scores
  message("computing  transcriptional heterogeneity scores")
  het <- compute.tITH(x = matrices)

  #--- adjust group names if subsampling was performed
  if (subsample>1) {
    het$state <- sub("-[^_]+$", "", het$state)
  }

  #--- return transcriptional heterogeneity scores
  return(het)


}

