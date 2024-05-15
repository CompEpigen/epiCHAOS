
#' @importFrom magrittr %>%
magrittr::`%>%`

#' @importFrom ggplot2 aes labs geom_bar theme_classic theme element_text
ggplot2::aes
ggplot2::labs
ggplot2::geom_bar
ggplot2::theme_classic
ggplot2::theme
ggplot2::element_text()

#' @importFrom stats lm reorder residuals
stats::lm
stats::reorder
stats::residuals

#' @importFrom jaccard jaccard
jaccard::jaccard

#' @importFrom stringr str_split
stringr::str_split

#' Compute epiCHAOS scores on a list of matrices
#'
#' @param x A list of binarised scATAC or other single cell matrices with equal rows.
#'
#' @return A dataframe containing epiCHAOS scores.
#' @export
#'
#' @examples
#' library(magrittr)
#' m1 <- sample(c(0,1), 5*10, prob=c(0.8,0.2), replace = TRUE) %>% matrix(5,10)
#' m2 <- sample(c(0,1), 5*10, prob=c(0.8,0.2), replace = TRUE) %>% matrix(5,10)
#' rownames(m1) <- rownames(m2) <- 1:5
#' colnames(m1) <- colnames(m2) <- paste0("cell", 1:10)
#' x <- list(m1=m1, m2=m2)
#' heterogeneity <- compute.eITH(x)
compute.eITH <- function(x) {

  # "het" will hold the heterogeneity scores for each condition
  het <- data.frame()
  for (d in names(x)) {

    message(paste0("calculating pairwise distances for ", d))

    #--- coerce to a dataframe
    temp <- x[[d]] %>% as.matrix() %>% as.data.frame()

    #--- move to the next calculation if an empty matrix is encountered
    if (sum(temp)==0) { next }

    #--- exclude any cells with zero counts
    temp <- temp[,colSums(temp)>0]
    x[[d]] <- temp

    #--- compute pairwise count-centred jaccard indices
    jac <- corrr::colpair_map(temp, jaccard::jaccard, center=T)
    jac <- apply(jac, 1, as.numeric)

    #--- commpute the mean of all pairwise differences
    het[d,"het"] <- mean(jac, na.rm = T)

  }

  message("Compiling epiCHAOS scores")

  #--- create dataframe to hold heterogeneity scores
  het$state <- rownames(het)

  #--- fit a linear regression model of raw epiCHAOS scores against the total count of the matrices snd take the residuals of the model as a count corrected score
  avg.count <- lapply(x, colMeans) %>% lapply(mean) %>% unlist()
  fit <- lm(het$het~avg.count)
  het$het.adj <- residuals(fit)

  #--- convert values to a range of 0,1
  het$het.raw <- (het$het - min(het$het)) / (max(het$het) - min(het$het))
  het$het.adj <- (het$het.adj - min(het$het.adj)) / (max(het$het.adj) - min(het$het.adj))

  #--- negate values so that higher score reflects higher heterogeneity
  het$het.raw <- 1-het$het.raw
  het$het.adj <- 1-het$het.adj


  #--- return a dataframe with raw and count-adjusted epiCHAOS scores for each group
  return(het[,c("het.adj", "het.raw", "state")])

}


#' Create group matrices
#'
#' @param counts A counts matrix representing e.g. a peaks-by-cells matrix in the case of scATAC data. The counts matrix will be binarised if not already.
#' @param meta Metadata including the column on which the data should be grouped, e.g. cluster, cell type.
#' @param colname The name of the column in the metadata on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in the provided metadata.
#' @param n The number of cells to subset for each group/cluster. This defaults to 100 cells.
#' @param m The minimum number of cells per group/cluster. This defaults to 20 cells. If fewer than m cells are found in a group/cluster, an epiCHAOS score is not computed for that group.
#' @param index The rows on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a vector of numerical indices or a vector or names corresponding to the rownames of interest in "counts".
#' @param binarise whether to binarise the matrices. This default to TRUE and is required for epiCHAOS analysis. Binarise can be set to false for heterogeneity analysis on scRNA-seq data for example.
#'
#' @return A list of binarised peaks/regions-by-cells matrices, one for each cluster/group of interest.
#' @export
#'
#' @examples
#' library(magrittr)
#' counts <- sample(c(0,1), 10*100, prob=c(0.8,0.2), replace = TRUE) %>% matrix(10,100)
#' rownames(counts) <- 1:nrow(counts)
#' colnames(counts) <- paste0("cell", 1:ncol(counts))
#' metadata <- data.frame(row.names=colnames(counts),
#' cluster=sample(c("group1", "group2", "group3"), ncol(counts), replace=TRUE))
#' matrices <- create.group.matrices(counts, metadata, "cluster", m=10, n=10)
create.group.matrices <- function(counts, meta, colname, n=100, m=20, index=NULL, binarise=TRUE) {

  meta$group <- meta[,colname]

  #--- if row indices are provided, subset the counts matrix for the specified rows
  if (is.null(index)) {index <- rownames(counts) }
  counts <- counts[index,]

  #--- create a list to hold counts matrices for each group/cluster
  matrices <- list()

  for (group in unique(meta$group)) {
    ids <- meta[meta$group==group, ] %>% rownames()

    #--- if the number of cells is smaller than a specified minimum, skip to the next group
    if (length(ids)<m) { next }

    matrices[[paste0("group-",group)]] <- counts[,ids] %>% as.matrix()

    #--- if more cells than selected n (defaults to 100 cells), downsample for n cells for that group
    if (length(ids)>n) { matrices[[paste0("group-",group)]] <- counts[,sample(ids, n)] %>% as.matrix() }

    #--- binarise counts if not already
    if (binarise == TRUE) {
      matrices[[paste0("group-",group)]][matrices[[paste0("group-",group)]]>1] <- 1
    }
  }

  #--- return list of matrices for epiCHAOS calculation
  return(matrices)
}


#' Compute epiCHAOS scores with per-chromosome count correction, for application to cancer datasets where large copy number alterations may result in differences in total counts per chromosome
#'
#' @param x A list of binarised single-cell matrices on which to compute epiCHAOS scores with per-chromosome counts adjustment. The rownames of each matrix should indicate the chromosome, start and end coordinates separated by ":", "-" or "_".
#'
#' @return A dataframe containing epiCHAOS scores.
#' @export
#'
#' @examples
#' library(magrittr)
#' m1 <- sample(c(0,1), 10*10, prob=c(0.8,0.2), replace = TRUE) %>% matrix(10,10)
#' m2 <- sample(c(0,1), 10*10, prob=c(0.8,0.2), replace = TRUE) %>% matrix(10,10)
#' rownames <- paste0(c("chr1", "chr2"), ":", seq(1000,100000, 10000), "-",
#' seq(2000,101000, 10000))
#' rownames(m1) <- rownames(m2) <- rownames
#' colnames(m1) <- colnames(m2) <- paste0("cell", 1:10)
#' x <- list(m1=m1, m2=m2)
#' heterogeneity <- compute.eITH.cancer(x)
compute.eITH.cancer <- function(x) {

  chromosomes <- rownames(x[[1]]) %>% str_split("-|_|:") %>% lapply("[", 1) %>% unlist() %>% unique()

  # "dists" will hold jaccard distances, "het.chr" will hold per-chromoosme heterogeneity scores, coverage will hold per chromosome counts
  het.chr <- dists <- coverage <- list()
  for (chr in paste0(chromosomes, ":")) {


    for (d in names(x)) {

      message(paste0("calculating pairwise distances for ", d, " ", chr))
      temp <- x[[d]] %>% as.matrix() %>% as.data.frame()

      #--- subset single cell matrix for selected chromosome
      temp <- temp[grepl(chr, rownames(temp)),]

      #--- compute pairwise Jaccard distances
      jac <- corrr::colpair_map(temp, jaccard, center=T)
      jac <- apply(jac, 1, as.numeric)

      # get mean of pairwise jaccard indices
      dists[[d]] <- colMeans(jac, na.rm = T) %>% mean()
      coverage[[d]] <- colMeans(temp) %>% mean()

    }

    #--- regress out counts
    fit <- lm(unlist(dists)~unlist(coverage))


    #--- create dataframe to hold heterogeneity scores
    het.chr[[chr]] <- (fit$residuals - min(fit$residuals)) / (max(fit$residuals) - min(fit$residuals))

  }

  het <-  plyr::rbind.fill(lapply(het.chr,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})) %>% colMeans(na.rm = T)
  het <- data.frame(het=het, state=names(het))

  #--- convert values to a range of 0,1
  het$het <- (het$het - min(het$het)) / (max(het$het) - min(het$het))  # Scale to a range of 0-1
  het$het <- 1-het$het


  return(het)

}


#' Compute epiCHAOS scores given a counts matrix and metadata
#'
#' @param counts A counts matrix representing e.g. a peaks-bycells matrix in the case of scATAC data. The counts matrix will be binarised if not already.
#' @param meta A dataframe containing metadata including at least the column on which the data should be grouped, e.g. cluster or cell type.
#' @param colname The name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".
#' @param n The number of cells to subset for each group/cluster. This defaults to 100.
#' @param index The rows on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a vector of numerical indices or a vector or names corresponding to the rownames of interest in "counts".
#' @param plot Defaults to FALSE. If TRUE, a bar plot will be returned as well as the epiCHAOS scores.
#' @param cancer Defaults to FALSE. If TRUE, the count correction step will be performed per-chromosome in order to take into account differences in coverage arising from large copy number alterations
#' @param subsample The number of times cells should be subsampled as "pseudoreplicates" for between-group comparisons. This defaults to 1. If > 1, "n" cells are subsampled from each group the specified number of times.
#'
#' @return A dataframe containing epiCHAOS scores for each group.
#' @export
#'
#' @examples
#' library(magrittr)
#' counts <- sample(c(0,1), 10*50, prob=c(0.8,0.2), replace = TRUE) %>% matrix(10,50)
#' rownames(counts) <- 1:nrow(counts)
#' colnames(counts) <- paste0("cell", 1:ncol(counts))
#' metadata <- data.frame(row.names=colnames(counts),
#' cluster=sample(c("group1", "group2", "group3"), ncol(counts), replace=TRUE))
#' heterogeneity <- epiCHAOS(counts, metadata, "cluster", subsample=10)
epiCHAOS <- function(counts, meta, colname=colnames(meta)[1], n=100, index=NULL, plot=F, cancer=F, subsample=1) {

  # bind variables
  state <- het.adj <- NULL

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

  message("computing epiCHAOS scores")

  #--- compute epiCHAOS scores
  if (cancer==T) {
    het <- compute.eITH.cancer(x = matrices)
  } else {
    het <- compute.eITH(x = matrices)
  }

  #--- adjust group names if subsampling was performed
  if (subsample>1) {
    het$state <- sub("-[^_]+$", "", het$state)
  }

  #--- if plot
  if (plot) {

    #--- return a barplot of epiCHAOS scores
    p <- ggplot2::ggplot(het, aes(x = reorder(state, het.adj), y = het.adj+0.01)) +
      geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.6)+
      labs(x="", y="epiCHAOS")+
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(list(het=het, barplot=p))

  } else {

    #--- return epiCHAOS scores
    return(het)
  }

}
