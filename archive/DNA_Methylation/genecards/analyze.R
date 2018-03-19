library(dplyr)
setwd("C:\\Git\\UGA\\DNA_Methylation\\genecards")
filenames <- list.files(".\\results", full.names = T)
gene_symbol <- list.files(".\\results")
gene_symbol <- sub(".txt$", "", gene_symbol)
results <- data.frame()
for (i in seq_along(filenames)) {
    df <- read.table(filenames[i], header = T, sep = "\t")
    df$gene_symbol <- gene_symbol[i]
    results <- rbind.data.frame(results, df)
}

results <- cbind(results$gene_symbol, results[ ,1:ncol(results)-1])
results$is_elite_enhancer <- ifelse(results$is_elite_enhancer == "True", 1, 0)
results <- results %>% arrange(gene_symbol, -is_elite_enhancer, -Enhancer_Score)
results$is_elite_enhancer <- ifelse(results$is_elite_enhancer == 1, "True", "False")
save(results, file = "./enhancer.RData")