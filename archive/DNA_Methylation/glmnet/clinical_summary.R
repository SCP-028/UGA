setwd("/lustre1/yz73026/array/clinical")
source("./GLM_func.R")
load("../annotation/annotation.RData")
tag.GLM <- "./GLM/"
tag.summary <- "./summary/"
GLM.normal <- list.files(tag.GLM, pattern = "normal")
projects <- sub("(^\\w+)\\..*", "\\1", GLM.normal)
summary.normal <- paste0(tag.summary, projects,
                         ".summary.normal.stage.RData")


getSummary <- function(glm) {
#   glm is a nested list of GLMs, raw results from clinical.R
    resultdf <- vector("list", length = length(glm))
    for (i in seq_along(glm)) {
        if (length(glm[[i]]) == 0 | is.null(glm[[i]])) {
            next
        }
        else {
            df <- buildDf(row.num=length(glm[[i]]), col.num=4)
            colnames(df) <- c("gene_symbol", "transcription_factor",
                              "lambda_1se", "lambda_min")
            df$gene_symbol <- sub("\\|.*$", "", names(glm[[i]]))
            df$transcription_factor <- lapply(glm[[i]], function(x) paste(x$TF, collapse = ";"))
            df$lambda_1se <- extractLambda(glm[[i]], func="lambda.1se")
            df$lambda_min <- extractLambda(glm[[i]], func="lambda.min")
            for (j in 2:5) {
                temp <- extractCoef(glm[[i]], j)
                df <- cbind.data.frame(df, temp)
            }
        }
        resultdf[[i]] <- df
    }
    names(resultdf) <- names(glm)
    return(resultdf)
}


detailSummary <- function(glm, predictor, resultdf, annot) {
    result <- vector("list", length=length(glm))
    for (i in seq_along(glm)) {
        if (length(glm[[i]]) == 0 | is.null(glm[[i]])) {
            next
        }
        else {
            df <- buildDf(row.num=length(glm[[i]]), col.num=4)
            colnames(df) <- c("gene_symbol",
                          paste0("coef", predictor, "_CpG"),
                          paste0("coef", predictor, "_%dev"),
                          paste0("coef", predictor, "_lambda"))
            rownames(df) <- names(glm[[i]])
            df$gene_symbol <- sub("\\|.*$", "", rownames(df))
            df[ ,3] <- resultdf[[i]][ ,paste0("coef", predictor, "_%dev")]
            df[ ,4] <- resultdf[[i]][ ,paste0("coef", predictor, "_lambda")]
            df[ ,2] <- getCpG(glm[[i]], predictor, resultdf[[i]])
            df <- calcGroup(df, resultdf[[i]], predictor, annot)
            df$Row.names <- NULL
            result[[i]] <- df
        }
        names(result) <- paste0("stage_", 1:4)
    }
    return(result)
}

groupCount <- function(coefn) {
    result <- matrix(nrow = 4, ncol = 7,
                     dimnames = list(paste0("stage_", 1:4),
                     c("1stExon", "3'UTR", "5'UTR",
                       "Body", "TSS1500", "TSS200",
                       "TFmethy")))
    for (i in seq_along(coefn)) {
        if (is.null(coefn[[i]])) {
            next
        }
        else {
            x <- colSums(coefn[[i]][ ,5:11 ])
            result[i, ] <- x
        }
    }
    result <- as.data.frame(result)
    result$stage <- rownames(result)
    return(result)
}

groupPlot <- function(groupCount_result_list, filename) {
#   Save a plot containing all four coefs.
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(reshape2))
    result <- data.frame()
    coef_name <- paste0("coef", 2:5)
    for (i in seq_along(groupCount_result_list)) {
        df <- melt(groupCount_result_list[[i]], id="stage")
        df$coef <- coef_name[i]
        result <- rbind.data.frame(result, df)
    }
    colnames(result) <- c("stage", "group", "value", "coef")
    ggplot(result, aes(x=group, y=value, colour=stage))+
    geom_path(aes(group=stage))+
    facet_wrap(~coef, scales = "free")+
    scale_color_brewer(palette = "RdYlBu")+
    ggsave(filename, dpi = 200, width = 20, height = 12, units = "in")
}


suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
registerDoMC(12)
foreach (i = seq_along(projects)) %dopar% {
    load(paste0(tag.GLM, GLM.normal[i]))
    results <- getSummary(normalGLM)  # List of 4 dfs
    coef2 <- detailSummary(normalGLM, 2, results, annot)
    coef3 <- detailSummary(normalGLM, 3, results, annot)
    coef4 <- detailSummary(normalGLM, 4, results, annot)
    coef5 <- detailSummary(normalGLM, 5, results, annot)

    save(results, coef2, coef3, coef4, coef5, file = summary.normal[i])
    groupCount_result_list <- list(coef2, coef3, coef4, coef5)
    groupCount_result_list <- lapply(groupCount_result_list, function(x) groupCount(x))
    groupPlot(groupCount_result_list, filename=paste0("./plot/", projects[i], "_normal.png"))
    print(paste0("Finished initial summary of ", projects[i]))
    rm(coef2, coef3, coef4, coef5)
}
