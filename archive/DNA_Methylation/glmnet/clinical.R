setwd("/lustre1/yz73026/array")
source("./clinical/GLM_func.R")

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
load("./annotation/sample2clinical.RData")

library(foreach, quietly = T, warn.conflicts = F)
library(doMC, quietly = T, warn.conflicts = F)
registerDoMC(12)

stages <- c("^Stage I[AB]?$", "^Stage II[ABC]?$", "Stage III", "Stage IV")


fitGLM <- function(methy, isoform, annotation, trans){
    library(glmnet, quietly = T, warn.conflicts = F)
    library(dplyr, quietly = T, warn.conflicts = F)
    library(data.table, quietly = T, warn.conflicts = F)
    annotation <- annotation %>% mutate_if(is.factor, as.character)
    trans <- trans %>% mutate_if(is.factor, as.character)
    annotation <- as.data.table(annotation)
    GLMdata <- vector("list", ifelse(!is.null(nrow(isoform)), nrow(isoform), 0))
    gene_name <- sub("\\|.*", "", rownames(isoform))
    tf_list <- lapply(gene_name, function(x) trans$TF[trans$Gene == x])
    CpG_list <- Map(c, gene_name, tf_list)
    CpG_list <- lapply(CpG_list, function(x) annotation[Gene_Name %in% unlist(x, use.names = F), IlmnID])
    for(i in seq_len(ifelse(!is.null(nrow(isoform)), nrow(isoform), 0))) {
        y <- t(isoform[i, , drop = F])
        x <- t(methy[rownames(methy) %in% CpG_list[[i]], , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1){
            next
        }
        else {
            tryCatch(GLMdata[[i]] <- list("Gene" = gene_name[i],
                       "TF" = tf_list[[i]],
                       "CpG site" = CpG_list[[i]],
                       "GLM" = suppressWarnings(cv.glmnet(x,y))),
                    error = function(e) e,
                    finally = next
            )
        }
    }
    names(GLMdata) <- rownames(isoform)
    GLMdata <- Filter(Negate(function(x) is.null(unlist(x))),
                      GLMdata)
    print(paste0("Made ", length(GLMdata), " GLMs."))
    return(GLMdata)
}


foreach (i = seq_along(projects)) %dopar% {
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
                if (!0 %in% dim(datan_stg) & !0 %in% dim(methyn_stg)) {
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
