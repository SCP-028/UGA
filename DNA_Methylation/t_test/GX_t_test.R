try(setwd("C:/Users/jzhou/Desktop"))
try(setwd("/lustre1/yz73026/array"))
filenames <- list.files("./mergeData/")
projects <- sub("\\.RData", "", filenames)
OUTPUT <- paste0("./t_test_expression/", projects, ".RData")
filenames <- list.files("./mergeData", full.names = T)

library(foreach)
library(doMC)
registerDoMC(12)
foreach (i = seq_along(projects)) %dopar% {
    load(filenames[i])
    df <- data.frame(gene_symbol=character(nrow(datan)),
                     p_value=numeric(nrow(datan)),
                     tumor_exp=factor(nrow(datan),
                                      levels = c("down-regulated",
                                                 "up-regulated",
                                                 "insignificant")))
    rm(methyn, methyt)
    df$gene_symbol <- rownames(datan)
    for (j in seq(nrow(datan))) {
        result <- t.test(datan[j, ], datat[j, ])
        names(result$estimate) <- NULL
        df$p_value[j] <- result$p.value
        if (result$p.value <= 0.05) {
            if (result$estimate[1] > result$estimate[2]) {
                df$tumor_exp[j] <- "down-regulated"
            }
            else {
                df$tumor_exp[j] <- "up-regulated"
            }
        }
        else {
            df$tumor_exp[j] <- "insignificant"
        }
        
    }
    save(df, file = OUTPUT[i])
    print(paste0("Finished project: ", projects[i]))
}