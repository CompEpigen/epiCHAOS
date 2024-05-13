## code to prepare `hemato_rna` dataset goes here
require(Seurat)

#--- rna
rna <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scRNA-Healthy-Hematopoiesis-191120.rds")
rna <- rna[,grepl("BMMC", rna$Group)]

# subset cell types
include <- c("01_HSC", "17_B", "19_CD8.N")
rna <- rna[,rna$BioClassification %in% include]

# subset 100 cells per cell-type
set.seed(10)
keep.ids <- c()
for (i in unique(rna$BioClassification)) {
  ids <- rna@colData[rna@colData$BioClassification==i,] %>% rownames() %>% sample(100)
  keep.ids <- c(keep.ids, ids)
}

rna <- rna[,keep.ids]


#--- create a seurat object, then normalise and scale the data
gex <- CreateSeuratObject(counts = rna@assays@.xData$data$counts, meta.data = data.frame(rna@colData))

# subset to 5000 genes to use only for example purposes
gex <- gex[sample(nrow(gex), 5000),]

#--- normalise and scale the data
gex <- NormalizeData(gex)
gex <- ScaleData(gex, features = rownames(gex))

hemato_rna <- list(data.scale=gex@assays$RNA$scale.data, meta=gex@meta.data)

# add data to epiCHAOS package
usethis::use_data(hemato_rna, overwrite = T)
