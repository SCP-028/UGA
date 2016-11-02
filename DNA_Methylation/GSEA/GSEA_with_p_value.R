# Uses glmnet results and perform a Gene Set Enrichment Analysis
# with Enrichment Scores and p-values calculated.

## Get glmnet results

## fit linear model and get R^2

## Resample using permutation test (shuffle order of gene expression)

## Run GSEA with the original model and the new 1000 models, 
## getting a p-value and a Enrichment Score

## Path and libraries
try(setwd("/lustre1/yz73026/array"))
rm(list = ls())
library(foreach)
library(doMC)
registerDoMC(12)

## Source codes and load data
source("./get_R2_of_glmnet.R")  # Code for messing with glmnet results

glmnet.normal <- list.files("./glmnet/", pattern = "\\.normal")
glmnet.tumor <- list.files("./glmnet/", pattern = "\\.tumor")
projects <- sub("(.*)\\.RData", "\\1", glmnet.normal)
coefNum <- 2:5
for (i in seq_along(coefNum)){
    
}
