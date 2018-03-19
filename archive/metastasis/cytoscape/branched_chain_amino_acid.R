library(dplyr)
setwd("~/data/BCAA")
load("../Cancers_Express.RData")
annot <- data.frame(
  `shared name`=c(
    "BCKA_ext (SLC16A1) BCKA_cyt",
    "BCKA_ext (SLC16A2) BCKA_cyt",
    "BCKA_ext (SLC16A4) BCKA_cyt",
    "BCKA_ext (SLC16A5) BCKA_cyt",
    "BCKA_ext (SLC16A6) BCKA_cyt",
    "BCKA_ext (SLC16A7) BCKA_cyt",
    "BCKA_cyt (BCAT1) BCAA_cyt",
    "BCAA_cyt (SLC7A5) BCAA_ext",
    "BCAA_cyt (SLC7A6) BCAA_ext",
    "BCAA_cyt (SLC7A7) BCAA_ext",
    "BCAA_cyt (SLC7A8) BCAA_ext",
    "BCAA_cyt (SLC43A1) BCAA_ext",
    "BCAA_cyt (SLC43A2) BCAA_ext",
    "BCAA_cyt (SLC43A3) BCAA_ext",
    "BCAA_cyt (?) BCAA_mc",
    "BCAA_mc (BCAT2) BCKA_mc",
    "BCKA_mc (BCKDHA) TCA",
    "BCKA_mc (BCKDHB) TCA",
    "BCKA_mc (DBT) TCA"
  ),
  source=c(
    "BCKA_ext", "BCKA_ext", "BCKA_ext", "BCKA_ext", "BCKA_ext", "BCKA_ext",
    "BCKA_cyt", "BCAA_cyt", "BCAA_cyt", "BCAA_cyt", "BCAA_cyt", "BCAA_cyt",
    "BCAA_cyt", "BCAA_cyt", "BCAA_cyt", "BCAA_mc", "BCKA_mc", "BCKA_mc","BCKA_mc"
  ),
  interaction=c(
    "SLC16A1", "SLC16A2", "SLC16A4", "SLC16A5", "SLC16A6", "SLC16A7",
    "BCAT1", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC43A1", "SLC43A2",
    "SLC43A3", "?", "BCAT2", "BCKDHA", "BCKDHB", "DBT"),
  target=c(
    "BCKA_cyt", "BCKA_cyt", "BCKA_cyt", "BCKA_cyt", "BCKA_cyt", "BCKA_cyt",
    "BCAA_cyt", "BCAA_ext", "BCAA_ext", "BCAA_ext", "BCAA_ext", "BCAA_ext",
    "BCAA_ext", "BCAA_ext", "BCAA_mc", "BCKA_mc", "TCA", "TCA",
    "TCA"
  ),
  check.names=F  # prevent changing `shared name` to `shared.name`
)

genelist <- c(
  "SLC16A1", "SLC16A2", "SLC16A4", "SLC16A5", "SLC16A6", "SLC16A7",
  "BCAT1", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC43A1", "SLC43A2",
  "SLC43A3", "BCAT2", "BCKDHA", "BCKDHB", "DBT")
genelist <- paste0(sub("(.*)", "(\\1$)", genelist), collapse="|")
result <- vector("list", length(All.DE))
for(i in seq_along(All.DE)) {
  df <- as.data.frame(All.DE[[i]])
  df <- df[grep(genelist, rownames(df)), ]
  df$interaction <- sub("^.*\\|(.*)$", "\\1", rownames(df))
  df$colors <- ifelse(df$p.fdr <= 0.05, ifelse(df$fold.change >= 1, "red", "blue"), "black")
  result[[i]] <- df
}
names(result) <- names(All.DE)
for(i in seq_along(result)) {
  df <- result[[i]]
  df <- left_join(annot, df, by="interaction")
  df$p.fdr[is.na(df$p.fdr)] <- 0
  df$fold.change[is.na(df$fold.change)] <- 1
  df$colors[is.na(df$colors)] <- "black"
  write.csv(df, file = paste0(names(result)[i], ".csv"), quote = F, row.names = F)
}
