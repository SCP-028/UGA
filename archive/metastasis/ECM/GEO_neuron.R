library(GEOquery)

setwd("C:/Users/jzhou/Desktop/")
df <- read.csv("geo_neuron.csv")
gseList <- as.character(df$Series)

results <- data.frame(GSE=character(), GSM=character(),
                      title=character(), description=character(),
                      method=character())

# Series
for (i in seq_along(gseList)) {
    gse <- getGEO(gseList[i], GSEMatrix = F, getGPL = F)
    # gsmList <- gse@header$sample_id
    gsmList <- lapply(gse@gsms, function(x) c(x@header$title, x@header$description))
    temp <- matrix(nrow = length(gsmList), ncol = 5)
    colnames(temp) <- c("GSE", "GSM", "title", "description", "method")
    for (j in seq_along(gsmList)) {
        temp[j,1:2] <- c(gseList[i], names(gsmList)[j])
        tryCatch(
                 temp[j,3:5] <- c(gsmList[[j]][1], gsmList[[j]][2], gsmList[[j]][3]),
                 error = function(e) print("Missing values, manually search for these"),
                 finally = next
                 )
    }
    results <- rbind.data.frame(results, temp)
}

results <- data.frame(GDS=character(), GSE=character(), description=character())
# DataSets
for (i in seq_along(gdsList)) {
    gds <- getGEO(gdsList[i], GSEMatrix = F, getGPL = F)
    gseDF <- cbind.data.frame(gds@dataTable@columns$sample, gds@dataTable@columns$description)
    gseDF$GDS <- gdsList[i]
    colnames(gseDF) <- c("GSM", "description", "GDS")
    gseDF <- cbind.data.frame(gseDF$GDS, gseDF[ ,1:2])
    results <- rbind.data.frame(results, gseDF)
}