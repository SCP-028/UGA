library(GEOquery)


makedf <- function(gds) {
    df <- getGEO(gds)
    df.table <- df@dataTable@table
    df.table <- df.table[grep("(^EFN[AB]\\d)|(EPH[AB]\\d{1,2})", df.table$IDENTIFIER, ignore.case = T), ]
    return(df.table)
}


df1 <- makedf("GDS5437")
df2 <- makedf("GDS3151")
df3 <- makedf("GDS2514")
