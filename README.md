# UGA
Code used at UGA for research.

## DNA methylation
This is where I put most effort into in 2016. Inside there's a workflow starting from data preprocessing, and moving on to unpaired t-tests, fitting generalized linear models, and finally GSEA (Gene Set Enrichment Analysis).

The folder also includes bisulfite sequencing data's methylation percentage calculation code. That includes genes and LINE-1 elements.

## Pathology
This is just some exercising work, downloading data from TCGA and organizing them. Deprecated.

## Warburg
Exercise work on modeling. Deprecated.

## Iron
Cancer occurrence rate ~ iron concentration in organs.

## Metastasis
Side project on cancer metastasis, mainly focusing on cadherin, a protein that holds cells together, and sialic acid, which is negatively charged and enriched on the cell membrane, therefore pushing cells away from each other.

Also contains DEG_func.R, the scripts for cleaning data, performing DEG analysis, and plotting the results.

### Classifier
Simple models for classifying different stages of breast cancer. Implemented with Python Sci-kit Learn.

### ECM
Differential Expressed Gene (DEG) analysis of extracellular matrix components.

### Repulsion
Co-expression and DEG analysis of cadherin and sialic acid genes.

### TF
Currently only scripts for calling gene symbols given the chromosome locations.

## pH
Calculating the overall charge of proteins.

## TCGA
Downloading all gene expression data from the GDC-TCGA database.

## UCP
Deprecated. DEG analysis and clustering.