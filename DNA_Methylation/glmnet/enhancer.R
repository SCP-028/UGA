setwd("/lustre1/yz73026/array")
source("./enhancer/GLM_func.R")

## Filenames
methy_files <- list.files("./methylation450", pattern = ".RData$")
isoform_files <- list.files("./expression", pattern = "isoform")
projects <- getProject(methy_files, isoform_files)
methy_files <- paste0("./methylation450/", projects, 
                      ".TNMethylation450.RData")
isoform_files <- paste0("./expression/", projects, 
                        ".v2.isoforms.normalized.RData")
result_files <- paste0("./clinical/GLM/", projects)

## Annotation files
load("./annotation/annotation.RData")
load("./annotation/cpgEnhancer2Targets.RData")

## Parallel settings
library(foreach, quietly = T, warn.conflicts = F)
library(doMC, quietly = T, warn.conflicts = F)
registerDoMC(12)

fitGLM <- function(methy, isoform, annotation, trans, enhancer){
    library(glmnet, quietly = T, warn.conflicts = F)
    library(dplyr, quietly = T, warn.conflicts = F)
    annotation <- annotation %>% mutate_if(is.factor, as.character)  # Convert factor to character
    GLMdata <- vector("list", nrow(isoform))
    gene_name <- sub("\\|.*", "", rownames(isoform))
    tf_list <- lapply(gene_name, function(x) as.character(trans$TF[trans$Gene == x]))
    eh_list <- lapply(gene_name, function(x) )
    CpG_list <- Map(c, gene_name, tf_list)
    CpG_list <- lapply(CpG_list, function(x) annotation$IlmnID[annotation$Gene_Name %in% unlist(x, use.names = F)])
    for(i in 1:nrow(isoform)){
        y <- t(isoform[i, , drop = F])
        x <- t(methy[rownames(methy) %in% CpG_list[[i]], , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1){
            next
        }
        else {
            df <- list("Gene" = gene_name[i],
                       "TF" = tf_list[[i]],
                       "CpG site" = CpG_list[[i]],
                       "GLM" = cv.glmnet(x,y, parallel = T))
            GLMdata[[i]] <- df
        }
    }
    names(GLMdata) <- rownames(isoform)
    GLMdata <- Filter(Negate(function(x) is.null(unlist(x))),
                      GLMdata)  # Remove empty entries
    print(paste0("Made ", length(GLMdata), " GLMs."))
    return(GLMdata)
}
for (i in seq_along(projects)) {
    print(paste0("Working on project: ", projects[i]))
    load(methy_files[i])
    if (ncol(datan) != 0 & ncol(datat) != 0) {      
        methyn <- datan
        load(isoform_files[i])
        if (!is.null(datan) & !is.null(datat)) {
            timer <- proc.time()
            rm(datat)
            datan <- cleanIsoform(datan, annot)
            rlist[methyn, datan] <- methyIsoform(methyn, datan, annot)
            normalGLM <- vector("list", length = length(stages))
            for (j in seq_along(stages)) {
                cases <- names(sample2stage[grep(stages[j], sample2stage)])
                datan_stg <- datan[ ,colnames(datan) %in% cases]
                methyn_stg <- methyn[ ,colnames(methyn) %in% cases]
                if (sum(dim(datan_stg) == 0) + sum(dim(methyn_stg) == 0) == 0) {
                    normalGLM[[j]] <- fitGLM(methyn_stg, datan_stg, annot, tf)
                }
            }
            names(normalGLM) <- paste0("Stages ", 1:4)
            save(normalGLM, file = paste0(result_files[i], ".normalGLM.stage.RData"))
            print(proc.time() - timer)
            rm(datan, methyn, datan_stg, methyn_stg)
        }
        else {
            rm(datan, datat, methyn)
        }
    }
    else {
        rm(datan, datat)
    }
    print(paste0("Finished ", projects[i]))
    gc()
}
