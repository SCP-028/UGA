mergeResult <- function(dirpath) {
    filenames <- list.files(dirpath, pattern = "RData$")
    projects <- sub("\\.RData", "", filenames)
    filenames <- list.files(dirpath, pattern = "RData$", full.names = T)
    result <- data.frame()
    for (i in seq_along(projects)) {
        load(filenames[i])
        df$cancer_type <- projects[i]
        result <- rbind.data.frame(result, df)
    }
    return(result)
}

countOccur <- function(df, type="expression", save.plot=T) {
    library(dplyr, quietly = T, warn.conflicts = F)
    if (grepl("exp(ression)?", type)) {
        result <- df %>% group_by(tumor_exp, cancer_type) %>%
                         summarise(length(cancer_type)) %>%
                         arrange(cancer_type)
        colnames(result)[3] <- "count"
        if (save.plot) {
            library(ggplot2, quietly = T, warn.conflicts = F)
            ggplot(df, aes(tumor_exp))+
            geom_bar()+
            facet_wrap(~cancer_type, scales = "free")+
            ggtitle("Gene Expression")+
            ggsave(filename = "./expression.png", width = 16, height = 9)
        }
    }
    else if (grepl("methy(lation)?", type)) {
        result <- df %>% group_by(tumor_methy, cancer_type) %>%
                         summarise(length(cancer_type)) %>%
                         arrange(cancer_type)
        colnames(result)[3] <- "count"
        if (save.plot) {
            library(ggplot2, quietly = T, warn.conflicts = F)
            ggplot(df, aes(tumor_methy))+
            geom_bar()+
            facet_wrap(~cancer_type, scales = "free")+
            ggtitle("Gene Methylation")+
            ggsave(filename = "./methylation.png", width = 16, height = 9)
        }
    }
    else {
        print("Invalid type, please check.")
    }
    return(result)
}

compareGene <- function(predictor, df, tumor_type, coef) {
#   predictor: data frames in predictors.RData
#   df: data frame generated from a loop below
#   tumor_type: e.g. "BLCA"  "COAD"
#   coef: e.g. 2 "coef5"
    library(dplyr, warn.conflicts = F, quietly = T)
    if (grepl("coef", coef)) {
        next
    }
    else {
        coef <- paste0("coef", coef)
    }
    df <- df %>%
            filter(cancer_type == tumor_type) %>%
            arrange(gene_symbol)
    predictor <- predictor %>%
                          filter(cancer_type == tumor_type & coef_num == coef) %>%
                          arrange(gene_symbol)
    up_df <- df %>% filter(tumor_exp == "up-regulated")
    down_df <- df %>% filter(tumor_exp == "down-regulated")
    up_predictor <- predictor[predictor$gene_symbol %in% up_df$gene_symbol, ]
    down_predictor <- predictor[predictor$gene_symbol %in% down_df$gene_symbol, ]
    coef_num <- as.numeric(sub(".*(\\d)", "\\1", coef))
    result <- list(up = up_predictor, down = down_predictor,
                   up_num = nrow(up_predictor) / coef_num,
                   down_num = nrow(down_predictor) / coef_num)
    return(result)
}


## Merge t_test results into df
try(setwd("C:/Users/jzhou/Desktop/t_test"))
filenames <- list.files("./", pattern = "RData$")
projects <- sub("\\.RData", "", filenames)
library(dplyr, warn.conflicts = F, quietly = T)
result <- data.frame()

for (i in seq_along(filenames)) {
    load(filenames[i])
    df$cancer_type <- projects[i]
    df <- df %>% filter(tumor_exp != "insignificant")
    result <- rbind.data.frame(result, df)
}
df <- result %>% arrange(cancer_type, p_value)

## Compare df and glmnet results, see distribution of well-predicted genes
try(setwd("../"))
load("predictors.RData")
up_predictor <- data.frame()
down_predictor <- data.frame()
num_genes <- data.frame(cancer_type = factor(levels = projects),
                        coef_num = factor(levels = c("coef2", "coef3", "coef4", "coef5")),
                        predictor = factor(levels = c("all.normal", "all.tumor",
                                                      "def0.6.normal", "def0.6.tumor")),
                        up_gene_num = numeric(),
                        down_gene_num = numeric())
p_list <- list(all.normal = all.normal, all.tumor = all.tumor,
               def0.6.normal = def0.6.normal, def0.6.tumor = def0.6.tumor)
pb <- txtProgressBar(min = 1, max = 192, char = "#", style = 3)
for (i in seq_along(projects)) {
    for (j in 2:5) {
        for (k in seq_along(p_list)) {
            ijk <- 16 * i + 4 * j + k -24
            result <- compareGene(predictor = p_list[[k]], df = df,
                                  tumor_type = projects[i], coef = j)
            up_predictor <- rbind.data.frame(up_predictor, result$up)
            down_predictor <- rbind.data.frame(down_predictor, result$down)
            num_genes[ijk, ] <- c(projects[i], paste0("coef", j), names(p_list)[k],
                                  result$up_num, result$down_num)
            setTxtProgressBar(pb, ijk)
        }       
    }
}
num_genes$down_gene_num <- as.numeric(num_genes$down_gene_num)
num_genes$up_gene_num <- as.numeric(num_genes$up_gene_num)

plotGenes <- function(df, predictors) {
#   df example: num_genes[num_genes$coef_num == "coef5", ]
    library(dplyr, warn.conflicts = F, quietly = T)
    library(ggplot2, warn.conflicts = F, quietly = T)
    temp <- df %>% filter(predictor == predictors)
    ggplot(temp, aes(down_gene_num, up_gene_num))+
    geom_point(alpha = 0.8, aes(color = cancer_type))+
    geom_abline(intercept = 0, slope = 1)+
    facet_wrap(~coef_num, scales = "free")+
    # geom_smooth(method = "lm", se = F)+
    ggtitle(predictors)
    ggsave(filename = paste0("./t_test_plot/", predictors, ".png"))
}

# coeflist <- c("coef2", "coef3", "coef4", "coef5")
predictorlist <- c("all.normal", "all.tumor",
                   "def0.6.normal", "def0.6.tumor")
for (j in seq_along(predictorlist)) {
    plotGenes(num_genes, predictorlist[j])
}
save(up_predictor, down_predictor, num_genes, file = "./t_test_result.RData")