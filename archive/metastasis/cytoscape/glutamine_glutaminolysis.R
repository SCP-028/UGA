library(dplyr)
setwd("~/data/glutaminolysis")
load("../Cancers_Express.RData")
annot <- data.frame(
  `shared name`=c(
    "glutamine_ext (SLC1A5) glutamine",
    "glutamine_ext (SLC6A14) glutamine",
    "glutamine_ext (SLC6A19) glutamine",
    "glutamine_ext (SLC7A5) glutamine",
    "glutamine_ext (SLC7A6) glutamine",
    "glutamine_ext (SLC7A7) glutamine",
    "glutamine_ext (SLC7A8) glutamine",
    "glutamine_ext (SLC7A9) glutamine",
    "glutamine_ext (SLC38A1) glutamine",
    "glutamine_ext (SLC38A2) glutamine",
    "glutamine_ext (SLC38A3) glutamine",
    "glutamine_ext (SLC38A5) glutamine",
    "glutamine_ext (SLC38A7) glutamine",
    "glutamine_ext (SLC38A8) glutamine",
    "glutamine (GLS) glutamate",
    "glutamine (GLS2) glutamate",
    "glutamate (GLUD1) alpha-ketoglutarate",
    "glutamate (GLUD2) alpha-ketoglutarate",
    "glutamate (GPT) alpha-ketoglutarate",
    "glutamate (GOT1) alpha-ketoglutarate",
    "glutamate (GOT2) alpha-ketoglutarate",
    "alpha-ketoglutarate (OGDH) succinyl-CoA",
    "alpha-ketoglutarate (DLD) succinyl-CoA",
    "alpha-ketoglutarate (DLST) succinyl-CoA",
    "succinyl-CoA (SUCLG1) succinate",
    "succinyl-CoA (SUCLG2) succinate",
    "succinyl-CoA (SUCLA2) succinate",
    "succinate (SDHA) fumarate",
    "succinate (SDHB) fumarate",
    "succinate (SDHC) fumarate",
    "succinate (SDHD) fumarate",
    "fumarate (FH) malate",
    "malate (MDH1) oxaloacetate",
    "malate (MDH2) oxaloacetate",
    "oxaloacetate (CS) citrate",
    "malate (ME1) pyruvate",
    "malate (ME2) pyruvate",
    "malate (ME3) pyruvate",
    "citrate (ACLY) Ac-CoA",
    "Ac-CoA (FASN) fatty_acid_synthesis"
  ),
  source=c(
    "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine_ext",
    "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine_ext",
    "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine_ext", "glutamine",
    "glutamine", "glutamate", "glutamate", "glutamate", "glutamate", "glutamate",
    "alpha-ketoglutarate", "alpha-ketoglutarate", "alpha-ketoglutarate", "succinyl-CoA",
    "succinyl-CoA", "succinyl-CoA", "succinate", "succinate", "succinate", "succinate",
    "fumarate", "malate", "malate", "oxaloacetate", "malate", "malate", "malate",
    "citrate", "Ac-CoA"
  ),
  interaction=c(
    "SLC1A5", "SLC6A14", "SLC6A19", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC7A9",
    "SLC38A1", "SLC38A2", "SLC38A3", "SLC38A5", "SLC38A7", "SLC38A8", "GLS", "GLS2","GLUD1",
    "GLUD2", "GPT", "GOT1", "GOT2", "OGDH", "DLD", "DLST", "SUCLG1", "SUCLG2", "SUCLA2",
    "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "CS", "ME1", "ME2", "ME3", "ACLY",
    "FASN"
  ),
  target=c(
    "glutamine", "glutamine", "glutamine", "glutamine", "glutamine", "glutamine", "glutamine",
    "glutamine", "glutamine", "glutamine", "glutamine", "glutamine", "glutamine", "glutamine",
    "glutamate", "glutamate", "alpha-ketoglutarate", "alpha-ketoglutarate",
    "alpha-ketoglutarate", "alpha-ketoglutarate", "alpha-ketoglutarate", "succinyl-CoA",
    "succinyl-CoA", "succinyl-CoA", "succinate", "succinate", "succinate", "fumarate",
    "fumarate", "fumarate", "fumarate", "malate", "oxaloacetate", "oxaloacetate","citrate",
    "pyruvate", "pyruvate", "pyruvate","Ac-CoA", "fatty_acid_synthesis"
  ),
  check.names=F  # prevent changing `shared name` to `shared.name`
)

genelist <- c(
  "SLC1A5", "SLC6A14", "SLC6A19", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC7A9",
  "SLC38A1", "SLC38A2", "SLC38A3", "SLC38A5", "SLC38A7", "SLC38A8", "GLS", "GLS2","GLUD1",
  "GLUD2", "GPT", "GOT1", "GOT2", "OGDH", "DLD", "DLST", "SUCLG1", "SUCLG2", "SUCLA2",
  "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "CS", "ME1", "ME2", "ME3", "ACLY",
  "FASN"
)
genelist <- paste0(sub("(.*)", "(\\\\\\|\\1$)", genelist), collapse="|")
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
