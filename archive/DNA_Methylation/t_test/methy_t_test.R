try(setwd("C:/Users/jzhou/Desktop"))
try(setwd("/lustre1/yz73026/array"))
filenames <- list.files("./mergeData/")
projects <- sub("\\.RData", "", filenames)
OUTPUT <- paste0("./t_test_methy/", projects, ".RData")
filenames <- list.files("./mergeData", full.names = T)

library(foreach)
library(doMC)
registerDoMC(12)
foreach (i = seq_along(projects)) %dopar% {
    load(filenames[i])
    methyn <- methyn[rownames(methyn) %in% rownames(methyt), ]
    methyt <- methyt[rownames(methyt) %in% rownames(methyn), ]
    df <- data.frame(ilmnID=character(nrow(methyn)),
                     p_value=numeric(nrow(methyn)),
                     tumor_methy=factor(nrow(methyn),
                                      levels = c("hypomethylated",
                                                 "hypermethylated",
                                                 "insignificant")))
    rm(datan, datat)
    df$ilmnID <- rownames(methyn)
    for (j in seq(nrow(methyn))) {
        result <- t.test(methyn[j, ], methyt[j, ])
        names(result$estimate) <- NULL
        df$p_value[j] <- result$p.value
        if (result$p.value <= 0.05) {
            if (result$estimate[1] > result$estimate[2]) {
                df$tumor_methy[j] <- "hypomethylated"
            }
            else {
                df$tumor_methy[j] <- "hypermethylated"
            }
        }
        else {
            df$tumor_methy[j] <- "insignificant"
        }
        
    }
    save(df, file = OUTPUT[i])
    print(paste0("Finished project: ", projects[i]))
}