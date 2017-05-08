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
## annot, tf
load("./annotation/annotation.RData")
## sample2{survival, grade, meta, stage}
load("./annotation/sample2clinical.RData")

library(foreach, quietly = T, warn.conflicts = F)
library(doMC, quietly = T, warn.conflicts = F)
registerDoMC(19)

stages <- c("^Stage I[AB]?$", "^Stage II[ABC]?$", "Stage III", "Stage IV")
grades <- paste0("G", 1:4)

fitGLM <- function(methy, isoform, annotation, trans){
    library(glmnet, quietly = T, warn.conflicts = F)
    GLMdata <- vector("list", nrow(isoform))
    gene_name <- sub("\\|.*", "", rownames(y))
    for(i in 1:nrow(isoform)){
        # Gene expression
        y <- isoform[i, , drop = F]
        # Add transcription factor name(s)
        z <- as.character(trans$TF[trans$Gene == gene_name[i]])
        z <- c(gene_name[i], z)
        CpG <- as.character(annotation[annotation$Gene_Name %in% z,1])
        x <- t(methy[rownames(methy) %in% CpG, , drop = F])
        if (is.null(dim(x)) | ncol(x) <= 1){
            next
        }
        else {
            df <- list("Gene" = gene_name[i],
                       "TF" = z[z != gene_name[i]],
                       "CpG site" = CpG,
                       "GLM" = cv.glmnet(x,y, parallel = T))
            GLMdata[[i]] <- df
        }
    }
    names(GLMdata) <- rownames(isoform)
    GLMdata <- Filter(Negate(function(x) is.null(unlist(x))),
                      GLMdata)
    return(GLMdata)
}

foreach (i = seq_along(projects)) %dopar% {
    ## Load and clean methylation and isoform data
    print(paste0("Working on project: ", projects[i]))
    load(methy_files[i])
    if (ncol(datan) != 0 & ncol(datat) != 0) {
        timer <- proc.time()
        methyn <- datan
        methyt <- datat
        load(isoform_files[i])
        if (!is.null(datan) & !is.null(datat)) {
<<<<<<< HEAD
            datan <- cleanIsoform(datan, annot)
            datat <- cleanIsoform(datat, annot)
            ## Pair methylation and isoform data
            rlist[methyn, datan] <- methyIsoform(methyn, datan, annot)
            rlist[methyt, datat] <- methyIsoform(methyt, datat, annot)
=======
            datan <- cleanIsoform(datan)
            datat <- cleanIsoform(datat)
            ## Pair methylation and isoform data
            rlist[methyn, datan] <- methyIsoform(methyn, datan)
            rlist[methyt, datat] <- methyIsoform(methyt, datat)
>>>>>>> bd8b3532e754da2200e0b6d48631d1291bcdb873
            ## Separate cases into different stages / grades
            normalGLM <- list()
            for (j in seq_along(stages)) {
                cases <- names(sample2stage[grep(stages[j], sample2stage)])
                datan_stg <- datan[ ,colnames(datan) %in% cases]
                methyn_stg <- methyn[ ,colnames(methyn) %in% cases]
                normalGLM[[j]] <- fitGLM(methyn_stg, datan_stg, annot, tf)
            }
            names(normalGLM) <- paste0("stage_", 1:4)
            save(normalGLM, file = paste0(results[i], ".normalGLM.RData"))
        }
    }
}