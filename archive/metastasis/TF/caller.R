require(biomaRt)
try(setwd("C:/Users/jzhou/Desktop/TF/data"), silent=T)
try(setwd("/home/yizhou/data/TF/data"), silent=T)


getAnnot <- function(mart, location) {
    location <- split(location, seq(ceiling(length(location) / 5000)))
    df <- data.frame()
    for (i in seq_along(location)) {
        annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
                       filters = 'chromosomal_region',
                       values = location[i],
                       mart = mart)
        df <- rbind.data.frame(df, annot)
        message(paste("Finished part", i, "..."))
    }
    return(df)
}


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
df <- data.table::fread("ENCFF804QQN.bed", header = F, stringsAsFactors = F,
                        data.table = F)
colnames(df) <- c("chrom", "chromStart", "chromEnd", "name", "score",
                  "strand", "signalValue", "pValue", "qValue", "peak")
df <- df[df$qValue >= 2, ]
# + strand
chromLocation <- paste(sub('chr', '', df$chrom),
                       df$chromStart,
                       df$chromEnd,
                       1, sep = ":")

annotP <- getAnnot(mart, chromLocation)

# - strand
chromLocation <- paste(sub('chr', '', df$chrom),
                       df$chromStart,
                       df$chromEnd,
                       -1, sep = ":")

annotN <- getAnnot(mart, chromLocation)
annot <- rbind.data.frame(annotP, annotN)
