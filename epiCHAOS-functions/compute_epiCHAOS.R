
#--- required packages
#install.packages(c("corrr", "jaccard", "dplyr", "magrittr", "plyr", "ggplot2"))
library(magrittr)


#--- function to compute epiCHAOS on a list of binarised single cell epigenomics matrices. Input is a list of binarised scATAC datasets with equal rows.
compute.eITH <- function(x) {  
  
  # "het" will hold the heterogeneity scores for each condition
  het <- data.frame()
  for (d in names(x)) { 
    
    print(paste0("calculating pairwise distances for ", d))
    
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
  
  print("Compiling epiCHAOS scores")
  
  #--- create dataframe to hold heterogeneity scores
  het$state <- rownames(het)
  
  #--- fit a linear regression model of raw epiCHAOS scores against the total count of the matrices snd take the residuals of the model as a count corrected score
  avg.count <- lapply(x, colMeans) %>% lapply(mean) %>% unlist()
  fit <- lm(het$het~avg.count)
  het$mean.het <- residuals(fit)
  
  #--- convert values to a range of 0,1
  het$het.raw <- (het$het - min(het$het)) / (max(het$het) - min(het$het))
  het$mean.het <- (het$mean.het - min(het$mean.het)) / (max(het$mean.het) - min(het$mean.het))
  
  #--- negate values so that higher score reflects higher heterogeneity
  het$het.raw <- 1-het$het.raw
  het$mean.het <- 1-het$mean.het
  
  
  #--- return a dataframe with raw and count-adjusted epiCHAOS scores for each group
  return(het[,c("mean.het", "het.raw", "state")])
  
}



#- counts: a counts matrix representing e.g. a peaks by cells matrix in the case of scATAC data. The counts matrix will be binarised if not already.
#- meta: metadata including the column on which the data should be grouped, e.g. cluster or celltype.
#- colname: the name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".
#- n: the number of cells to subset for each group/cluster. This defults to 100.
#- m: the minimum number of cells per group/cluster. This defaults to 20. If fewer than m cells are found in a group/cluster, an epiCHAOS score is not computed.
#- index: the rows in counts on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a (i) vector of numerical indices, (ii) a vector or names corresponding to the rownames of interest in "counts"

create.group.matrices <- function(counts, meta, colname, n=100, m=20, index=NULL) {

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
    matrices[[paste0("group-",group)]][matrices[[paste0("group-",group)]]>1] <- 1
  }
  
  #--- return list of matrices for epiCHAOS calculation
  return(matrices)
}


#--- function to compute epiCHAOS with count correction per chromosome, for application to cancer datasets where large copy number alterations may result in differences in total counts per chromosome

#- input is a list of binarised matrices
#- rownames in each matrix are the chromosome, start and end coordinates separated by ":", "-" or "_"
compute.eITH.cancer <- function(x) {
  
  chromosomes <- rownames(x[[1]]) %>% str_split("-|_|:") %>% lapply("[", 1) %>% unlist() %>% unique()
  
  # "dists" will hold jaccard distances, "het.chr" will hold per-chromoosme heterogeneity scores, coverage will hold per chromosome counts
  het.chr <- dists <- coverage <- list()
  for (chr in paste0(chromosomes, ":")) {
    
    
    for (d in names(x)) { 
      
      print(paste0("calculating pairwise distances for ", d, " ", chr))
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
  
  require(plyr)
  het <-  rbind.fill(lapply(het.chr,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)})) %>% colMeans(na.rm = T)
  het <- data.frame(het=het, state=names(het))
  
  #--- convert values to a range of 0,1
  het$het <- (het$het - min(het$het)) / (max(het$het) - min(het$het))  # Scale to a range of 0-1
  het$het <- 1-het$het
  
  
  return(het)
  
}


#--- function to compute epiCHAOS if given a counts matrix and metadata

# - counts: a counts matrix representing e.g. a peaks by cells matrix in the case of scATAC data. The counts matrix will be binarised if not already.
# - meta: metadata including the column on which the data should be grouped, e.g. cluster or celltype.
# - colname: the name of the column in "meta" on which to group the data e.g. "cluster" or "celltype". If not specified this defaults the the first column in "meta".
# - n: the number of cells to subset for each group/cluster. This defults to 100.
# - index: the rows in counts on which to subset the counts matrix. If not provided the whole counts matrix will be used by default. Otherwise "index" may be specified as either a (i) vector of numerical indices, (ii) a vector or names corresponding to the rownames of interest in "counts"
# - plot: if true, a barplot will be returned as well as the epiCHAOS scores. Defaults to FALSE.
# - cancer: if true, the count correction step will be performed per-chromosome in order to take into account differences in coverage arising from large copy number alterations
# - subsample: if > 1, the specified number of subsamples of cells are selected from each group for computation. 

epiCHAOS <- function(counts, meta, colname=colnames(meta)[1], n=100, index=NULL, plot=F, cancer=F, subsample=1) {
  
  #--- create per-group matrices
  print("creating group matrices")
  
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
  
  print("computing epiCHAOS scores")
  
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
    p <- ggplot2::ggplot(het, aes(x = reorder(state, mean.het), y = mean.het+0.01)) +
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



