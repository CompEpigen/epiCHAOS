## code to prepare `encode_tfbs` dataset goes here

#--- the lola encode tfbs downloaded from
lola <- get(load("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/encode_tfbs.RData"))
index <- read.table("/omics/groups/OE0219/internal/KatherineK/Genomic_Annotations/LOLA/hg19/encode_tfbs/index.txt", header = T)
names(lola) <- index$filename

# subset to use as example
lola <- lola[1:10]
index <- index[index$filename %in% names(lola),]

encode_tfbs <- list(ranges=lola, anno=index)


# add data to epiCHAOS package
usethis::use_data(encode_tfbs, overwrite = T)
