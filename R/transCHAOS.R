
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
    print(paste0("computing pairwise distances for ", i))
    temp <-  x[[i]] %>% t() %>% dist() %>% as.matrix() %>% colMeans()

    #--- check for and remove outliers in dists[[i]]
    bound <- mean(temp) + 3*sd(temp)
    dists[[i]] <- temp[temp<bound]

  }

  #--- create dataframe to hold heterogeneity scores
  print("Compiling scores to dataframe")
  het <- data.frame(het=unlist(dists))  # this assumes that all datasets have an equal number of cells
  state <- c()
  for (cond in names(dists)) {
    state <- c(state, rep(cond, length(dists[[cond]])))
  }

  het$state <- state

  for (state in unique(het$state)) {
    het$mean.het[het$state==state] <- mean(het$het[het$state==state])
  }

  het <- het[,c("state", "mean.het")] %>% unique()
  het <- na.omit(het)

  #--- transform to a range of 0-1
  het$mean.het <- (het$mean.het - min(het$mean.het)) / (max(het$mean.het) - min(het$mean.het))

  #--- return heterogeneity score
  return(het)
}
