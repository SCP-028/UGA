setwd("C:/Users/yz73026/Desktop/feb2018_paper_boxplots/tryptophan")
library(jsonlite)
cpds <- fromJSON("./Tryptophan_reactions.json")

result <- vector("list", length(cpds))
for (i in seq_along(cpds)) {
    cpd <- cpds[[i]]
    enzymes <- lapply(cpd, function(x) x$enzyme)
    Filter(function(x) !grepl("\\.", x), enzymes)
    result[[i]] <- unique(Filter(function(x) !grepl("\\.", x), unlist(enzymes)))
}
names(result) <- names(cpds)
