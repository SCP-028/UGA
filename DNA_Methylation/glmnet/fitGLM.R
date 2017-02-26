# Pre-process methy and GX datasets of different cancers,
# and fit GLM model.
try(setwd("/home/yi/data/DNA_Methylation/Methy_array/RData/"))
try(setwd("/lustre1/yz73026/array"))
source("./glmnet_TF/GLM_func.R")  # Functions for manipulating data and fitting GLM

# Find files to load
methyFiles <- list.files("./methylation450", pattern = "\\.RData")
expressionFiles <- list.files("./expression", pattern = "\\.RData")
projects <- getProject(methyFiles, expressionFiles)
methyFiles <- paste0("./methylation450/", projects, 
                ".TNMethylation450.RData")
expressionFiles <- paste0("./expression/", projects, 
                ".update1.genes.normalized.RData")
resultFiles <- paste0("./glmnet_mTF/", projects)

# Remove duplicated IlmnID, Gene name, and Gene group
# source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")
list[annot, tf] <- readAnnot()

library(foreach)
library(doMC)
registerDoMC(length(projects) / 2)  # 14 / 2 = 7
foreach (i = seq_along(projects)) %dopar% {
#   Different tumor types
    print(paste0("Working on project: ", projects[i]))
    load(methyFiles[i])
    if(ncol(datan) != 0 && ncol(datat) != 0){
#       Some data seem to be missing
        timer <- proc.time()
        methyn <- datan
        methyt <- datat
        load(expressionFiles[i])
        datan <- convertName(datan)
        datat <- convertName(datat)
        list[methyn, datan] <- pairNames(methyn, datan, annot)
        list[methyt, datat] <- pairNames(methyt, datat, annot)
        normalGLM <- fitGLM(methyn, datan, annot, tf)
        save(normalGLM, file = paste0(resultFiles[i],
                                    ".normalGLM.RData"))
        tumorGLM <- fitGLM(methyt, datat, annot, tf)
        save(tumorGLM, file = paste0(resultFiles[i],
                                    ".tumorGLM.RData"))
        print(proc.time() - timer)
        print(paste0("Project ", projects[i], " finished!"))
    }
}
