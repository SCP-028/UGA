setwd("~/data/DNA_Methylation/LINE1/Data/GSE81227")
filepath <- list.files("./", pattern = "\\.(txt)$", recursive = T)
results <- data.frame("Sample" = NA, "mCG" = NA, "Methy_Percent" = NA)

meanMethy <- function(filepath)
{
    for (i in 1:length(filepath))
    {
        filename <- sub("\\_.*\\.txt", "", filepath[i])
        df <- read.table(filepath[i], header = T)
        mCG <- mean(as.numeric(df$num_meth_CpG))
        percent <- mean(as.numeric(df$percent_meth)) * 0.01
        temp <- c(filename, mCG, percent)
        results <- rbind(results, temp)
    }
    results <- results[-1, ]
    rownames(results) <- 1:nrow(results)
    return(results)
}

df <- meanMethy(filepath)
write.table(df, "../results/GSE81227.txt")