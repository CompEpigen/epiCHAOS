## code to prepare `hemato` dataset goes here

# if(require("GenomicRanges")) { library(GenomicRanges) }
# if(require("S4Vectors")) { library(GenomicRanges) }

# single cell atac-seq data from human hematopoiesis downloaded from Granja et al. 2019: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7258684/
# peaks-by-cells data downloaded from github as an rds file: https://github.com/GreenleafLab/MPAL-Single-Cell-2019; https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds
atac <- readRDS("/omics/groups/OE0219/internal/KatherineK/data/scATAC/Granja2019/scATAC-Healthy-Hematopoiesis-191120.rds")

# data subsetted for bone-marrow-derived cells
atac <- atac[,grepl("BMMC", atac@colData$Group)]

# subset random 10000 peaks
#set.seed(10)
atac <- atac[sample(nrow(atac), 10000),]

# subset celltypes to have only progenitor cells, monocytes, B, T and NK cells
atac@colData$BioClassification %>% table()
exclude <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso", "09_pDC", "10_cDC", "13_Unk",  "14_Unk", "18_Plasma", "26_Unk", "08_GMP.Neut")
atac <- atac[,!atac$BioClassification %in% exclude]
atac$BioClassification %>% table()


# extract counts matrix and metadata
counts <- atac@assays$data$counts
meta <- atac@colData
ranges <- atac@rowRanges
hemato <- list(counts=counts, meta=meta, ranges=ranges)

# add data to epiCHAOS package
usethis::use_data(hemato, overwrite = T)
