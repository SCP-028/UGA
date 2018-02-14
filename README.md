# UGA
Code used at UGA for research.

## DNA methylation
This is where I put most effort into in 2016, resulting in the paper ["Transcription regulation by DNA methylation under stressful conditions in human cancer"](https://link.springer.com/article/10.1007/s40484-017-0129-y).

The folder also includes code for calculating the methylation percentage of bisulfite sequencing data. That includes genes and LINE-1 elements.

## Glycosylation
Some duplicated code for plotting expression information in a specific pathway. The `pathview` part might come useful in the future.

## Histone acetylation
Some unfinished work on solving ODEs in a acetylation system. Might port to Python when I have enough time.

## Metastasis
Side project on cancer metastasis, mainly focusing on cadherin, a protein that holds cells together, and sialic acid, which is negatively charged and enriched on the cell membrane, therefore pushing cells away from each other.

### Classifier
Simple models for classifying different stages of breast cancer. Implemented with Python Sci-kit Learn.

### CytoScape
R code framework for generating csv files that can be imported into CytoScape as a network for visualization.

### ECM
Differential Expressed Gene (DEG) analysis of extracellular matrix components.

### Repulsion
Co-expression and DEG analysis of cadherin and sialic acid genes.

### TF
Currently only scripts for calling gene symbols given the chromosome locations.

### DEG_func.R
Scripts for cleaning data, performing DEG analysis, and plotting the results.

## Misc
A bit of this and a bit of that. Archived.

### Iron
Cancer occurrence rate ~ iron concentration in organs.

### Pathology
This is just some exercising work, downloading data from TCGA and organizing them. Deprecated.

### Scripts
Correlation test of genes.

### TCGA
Downloading all gene expression data from the GDC-TCGA database.

### UCP
Deprecated. DEG analysis and clustering.

### Warburg
Exercise work on modeling. Deprecated.

## Protein pKa
Calculating protein overall charge with [MCCE](http://www.sci.ccny.cuny.edu/~mcce/index.php). Note that this software only works on Intel CPUs. 
