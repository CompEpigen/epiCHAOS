
# EpiCHAOS

EpiCHAOS is a metric for calculating cell-to-cell heterogeneity in single-cell epigenomics data. EpiCHAOS can be applied to any single cell epigenomics dataset (e.g. scATAC-seq, scChIP-seq, scNMT-seq, scTAM-seq) in binarised matrix form.

EpiCHAOS scores are assigned per cluster (or otherwise defined group of single cells, e.g. cell type, treatment condition). Scores will range from 0-1 where higher scores indicate higher cell-to-cell heterogeneity. 

For details and examples of applications please refer to the epiCHAOS manuscript:
Kelly, K., Scherer, M., Braun, M.M. et al. EpiCHAOS: a metric to quantify epigenomic heterogeneity in single-cell data. Genome Biol 25, 305 (2024). https://doi.org/10.1186/s13059-024-03446-w

#### Description of epiCHAOS calculation
1. Single-cell matrix is binarised if not already
2. Pairwise chance-corrected Jaccard indices are calculated between all cells in each group/cluster
3. The mean of all pairwise indices is computed per cluster as a raw heterogeneity score
4. Raw heterogeneity scores are adjusted for differences in sparsity between groups by fitting a linear regression model
5. Scores are normalised to a range of 0-1 and negated so that a higher score indicates higher cell-to-cell heterogeneity
   <img width="778" alt="epiCHAOS_schematic" src="https://github.com/CompEpigen/epiCHAOS/assets/61455651/0fdc19e5-7b50-4475-98b0-4ece1f3762a0">

#### Install the epiCHAOS R package
```
library(devtools)
install_github("CompEpigen/epiCHAOS")
```
#### Calculation of epiCHAOS scores
Calculation of epiCHAOS scores requires (i) a single cell epigenomics dataset in binarised matrix form, e.g. a peaks-by-cells or tiles-by-cells matrix in the case of scATAC-seq data, (ii) cell metadata annotation matching column names to clusters, cell types or other groups on which heterogenetiy scores should be computed.

The main function for computing epiCHAOS scores is "epiCHAOS()" which will take as input a single-cell counts matrix together with cell annotation metadata specifying the group on which cell-to-cell heterogeneity should be compared.

```
# compute epiCHAOS scores
heterogeneity <- epiCHAOS(counts = counts, meta = metadata, colname = "seurat_clusters", n = 50, index = NULL, plot = F, cancer = F, subsample = 5)
```
#### Demonstration using the hematopoietic example dataset
EpiCHAOS includes a small example single-cell ATAC-seq dataset which we will use to demonstrate its functionality. The example dataset is derived from the [Granja et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7258684/) dataset from the human hematopoietic system.

```
# example scATAC-seq data
counts <- hemato_atac$counts
rownames(counts) <- 1:nrow(counts)
metadata <- hemato_atac$meta %>% data.frame()

# confirm that "counts" is a matrix with names rows, and whose colnames correspond to the rownames of "metadata"
head(counts)
colnames(counts) %in% rownames(metadata)

# compute epiCHAOS scores for each cell type (cell type is found in the "BioClassification" column of "metadata")
heterogeneity <- epiCHAOS(counts = counts, meta = metadata, colname = "BioClassification", n = 50, index = NULL, plot = F, cancer = F, subsample = 1)

# inspect epiCHAOS scores. Note that "het.adj" is the epiCHAOS score, which has been adjusted for sparsity. The "het.raw" colunm is also returned to give you an idea of how the correction for sparsity affected the resulting scores. 
heterogeneity

# plot heterogeneity scores
ggplot(heterogeneity, aes(x = reorder(state, het.adj), y = het.adj+0.01)) +
  geom_bar(stat="identity", position = "dodge", alpha=0.8, width = 0.6)+
  labs(x="", y="epiCHAOS")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```
