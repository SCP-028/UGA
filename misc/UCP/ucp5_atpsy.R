## Multiple regression of UCPs against ATP synthases.
setwd("C:/Users/jzhou/Desktop/SLC25A14")
library(glmnet)
source("https://raw.githubusercontent.com/ggrothendieck/gsubfn/master/R/list.R")

file_name <- list.files("./", pattern="RData$")
projects <- sub("(.*)\\.up.*$", "\\1", file_name)

results <- vector("list", length(projects) * 2)
names(results) <- as.vector(outer(projects, c("_normal", "_tumor"), paste0))
names(results) <- names(results)[order(names(results))]
res.param <- matrix(nrow = length(results), ncol = 2, dimnames = list(names(results), c("coef_num", "dev_ratio")))


ucp5.t.test <- function(datan, datat) {
    x <- datan[grep("SLC25A14", rownames(datan)), ]
    y <- datat[grep("SLC25A14", rownames(datat)), ]
    if (class(x) == "numeric") {  # Got hits, otherwise would be empty matrix
        z <- t.test(x, y, na.rm=T)
        return(z$p.value <= 0.01 & z$estimate[1] < z$estimate[2])
    }
    else {
        return(FALSE)
    }
}


fit.glm <- function(ucp, atp, alpha=0) {
    # x <- t(rbind(ucp, atp))[ ,-1, drop=F]
    # y <- t(ucp[1, , drop=F])
    x <- rbind(atp, ucp)
    y <- t(x[1, , drop=F])
    x <- t(x[-1, , drop=F])
    result <- cv.glmnet(x, y, alpha=alpha)
    dev_ratio <- result$glmnet.fit$dev.ratio[which(result$glmnet.fit$lambda == result$lambda.min)]
    result <- as.matrix(coef(result, s="lambda.min"))
    result <- result[order(abs(result), decreasing=T), , drop=F]
    colnames(result) <- colnames(y)
    return(list(result, ncol(x), dev_ratio))
}


for (i in seq_along(projects)) {
    load(file_name[i])
    # if (ucp5.t.test(datan, datat)) {
        x.normal <- datan[grep("(^UCP1\\|)|(^UCP2\\|)|(^UCP3\\|)|(^SLC25A27\\|)|(^SLC25A14\\|)", rownames(datan)), , drop=F]
        # y.normal <- datan[grep("(^ATPIF1\\|)|(^ATP5A1\\|)|(^ATP5B\\|)|(^ATP5C1\\|)|(^ATP5D\\|)|(^ATP5E\\|)|(^ATP5F1\\|)|(^ATP5G1\\|)|(^ATP5G2\\|)|(^ATP5G3\\|)|(^ATP5H\\|)|(^ATP5I\\|)|(^ATP5J\\|)|(^ATP5J2\\|)|(^ATP5L\\|)|(^ATP5O\\|)|(^MT-ATP6\\|)|(^MT-ATP8\\|)", rownames(datan)), , drop=F]
        y.normal <- datan[grep("(^ATP5A1\\|)|(^ATP5D\\|)", rownames(datan)), , drop=F]
        x.tumor <- datat[grep("(^UCP1\\|)|(^UCP2\\|)|(^UCP3\\|)|(^SLC25A27\\|)|(^SLC25A14\\|)", rownames(datat)), , drop=F]
        # y.tumor <- datat[grep("(^ATPIF1\\|)|(^ATP5A1\\|)|(^ATP5B\\|)|(^ATP5C1\\|)|(^ATP5D\\|)|(^ATP5E\\|)|(^ATP5F1\\|)|(^ATP5G1\\|)|(^ATP5G2\\|)|(^ATP5G3\\|)|(^ATP5H\\|)|(^ATP5I\\|)|(^ATP5J\\|)|(^ATP5J2\\|)|(^ATP5L\\|)|(^ATP5O\\|)|(^MT-ATP6\\|)|(^MT-ATP8\\|)", rownames(datat)), , drop=F]
        y.tumor <- datat[grep("(^ATP5A1\\|)|(^ATP5D\\|)", rownames(datat)), , drop=F]
        rownames(x.normal) <- sub("(.*)\\|.*$", "\\1", rownames(x.normal))
        rownames(x.tumor) <- sub("(.*)\\|.*$", "\\1", rownames(x.tumor))
        rownames(y.normal) <- sub("(.*)\\|.*$", "\\1", rownames(y.normal))
        rownames(y.tumor) <- sub("(.*)\\|.*$", "\\1", rownames(y.tumor))
        if (nrow(x.normal) > 1 & nrow(y.normal) != 0) {
            list[results[[2*i-1]], res.param[2*i-1, 1], res.param[2*i-1, 2]] <- fit.glm(x.normal, y.normal, alpha=1)            
        }
        if (nrow(x.tumor) > 1 & nrow(y.tumor) != 0) {
            list[results[[2*i]], res.param[2*i, 1], res.param[2*i, 2]] <- fit.glm(x.tumor, y.tumor, alpha=1)            
        }

        # results[[2*i-1]] <- dfn <- cor(t(transporter_normal))
        # results[[2*i]] <- dft <- cor(t(transporter_tumor))
        # png(filename=paste0("./plot/ATPsy/", names(results)[2*i-1], ".png"), width = 1024, height = 1024, units = "px")
 
        # dev.off()
    # }
}

df <- character()
for (i in seq_along(results)) {
    df <- union(df, rownames(results[[i]]))
}

df <- data.frame(row.names = df)
library(dplyr)
for (i in seq_along(results)) {
    df <- left_join(df, as.data.frame(results[[i]]), by = c(rownames(df), rownames(results[[i]])))
}
colnames(df) <- names(results)

save(results, res.param, file="C:/Users/jzhou/Desktop/ucp5_atp.RData")
