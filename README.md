
# EpiCHAOS

EpiCHAOS is a metric for calculating cell-to-cell heterogeneity in single-cell epigenomics data. EpiCHAOS can be applied to any single cell epigenomics dataset (e.g. scATAC-seq, scChIP-seq, scNMT-seq, scTAM-seq) in binarised matrix form.

EpiCHAOS scores are assigned per cluster (or otherwise defined group of single cells, e.g. cell type, treatment condition). Scores will range from 0-1 where higher scores indicate higher cell-to-cell heterogeneity. 

For details please refer to the epiCHAOS preprint: https://www.biorxiv.org/content/10.1101/2024.04.24.590899v1. "analysis-scripts" contains R scripts used for all analyses related to this publication.

#### Description of epiCHAOS calculation
1. Single-cell matrix is binarised if not already
2. Pairwise chance-corrected Jaccard indices are calculated between all cells in each group/cluster
3. The mean of all pairwise indices is computed per cluster as a raw heterogeneity score
4. Raw heterogeneity scores are adjusted for differences in total counts between groups by fitting a linear regression model and taking the residuals of the model as an adjusted score
5. Scores are normalised to a range of 0-1 and subtracted from one so that a higher score indicates higher cell-to-cell heterogeneity
   <img width="778" alt="epiCHAOS_schematic" src="https://github.com/CompEpigen/epiCHAOS/assets/61455651/0fdc19e5-7b50-4475-98b0-4ece1f3762a0">

#### install the epiCHAOS R package
´´´
library(devtools)
install_github("colomemaria/epiAneufinder")
´´´

#### Calculation of epiCHAOS scores
Calculation of epiCHAOS scores requires (i) a single cell epigenomics dataset in binarised matrix form, e.g. a peaks-by-cells or tiles-by-cells matrix in the case of scATAC-seq data, (ii) cell annotation matching column names to clusters, cell types or other groups on which heterogenetiy scores should be computed.

The main function for computing epiCHAOS scores is "epiCHAOS()" which takes a counts matrix and cell annotation metadata as input and returns a dataframe of epiCHAOS scores for each group/cluster. 

´´´
# compute epiCHAOS scores
heterogeneity <- epiCHAOS(counts = counts, meta = metadata, colname = "BioClassification", n = 50, index = NULL, plot = F, cancer = F, subsample = 1)
´´´
