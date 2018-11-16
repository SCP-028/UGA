setwd("C:/Users/jzhou/Desktop/SLC25A14")
file_name <- list.files("./", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)
library(glmnet)
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")

fit.glm <- function(syndf, condf) {
#   Fit GLM for 2 input df.
    x <- t(rbind(syndf, condf))[ ,2:(nrow(syndf)+nrow(condf)), drop=F]
    y <- t(syndf[1, , drop=F])
    result <- cv.glmnet(x, y, alpha=0)
    dev_ratio <- result$glmnet.fit$dev.ratio[which(result$glmnet.fit$lambda == result$lambda.min)]
    result <- as.matrix(coef(result, s="lambda.min"))
    colnames(result) <- colnames(y)
    return(list(result, ncol(x), dev_ratio))
}

result.normal <- vector("list", length(projects))
result.tumor <- vector("list", length(projects))
predictor.normal <- vector("numeric", length(projects))
predictor.tumor <- vector("numeric", length(projects))
dev.normal <- vector("numeric", length(projects))
dev.tumor <- vector("numeric", length(projects))


`
for (i in seq_along(projects)) {
    load(file_name[i])
    syn.normal <- datan[grep("(^PDHA1\\|)|(^PDHB\\|)|(^CPT2\\|)|(^CPT1A\\|)|(^CPT1B\\|)|(^ACSS1\\|)", rownames(datan)), , drop=F]
    syn.tumor <- datat[grep("(^PDHA1\\|)|(^PDHB\\|)|(^CPT2\\|)|(^CPT1A\\|)|(^CPT1B\\|)|(^ACSS1\\|)", rownames(datat)), , drop=F]
    con.normal <-  datan[grep("(^ACOT12\\|)|(\\|1431$)", rownames(datan)), , drop=F]
    con.tumor <- datat[grep("(^ACOT12\\|)|(\\|1431$)", rownames(datat)), , drop=F]
    list[result.normal[[i]], predictor.normal[i], dev.normal[i]] <- fit.glm(syn.normal, con.normal)  # Coef of genes from syn.normal should be negative
    list[result.tumor[[i]], predictor.tumor[i], dev.tumor[i]] <- fit.glm(syn.tumor, con.tumor)
}
names(predictor.normal) <- names(dev.normal) <- names(predictor.tumor) <- names(dev.tumor) <- names(result.normal) <- names(result.tumor) <- projects



save(result.normal, result.tumor, predictor.normal, predictor.tumor, dev.normal, dev.tumor, file = "C:/Users/jzhou/Desktop/acetylation_result.RData")