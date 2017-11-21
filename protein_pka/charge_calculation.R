# Calculate protein net charge at pH = 7, then we can know
# how many protons were contributed by the protein.
# Need http://www.uniprot.org/uploadlists/ to convert IDs.
setwd("~/data/protein_pka")
require(dplyr)
require(tidyr)
require(ggplot2)

df <- data.table::fread("./pKa.csv", stringsAsFactors = F,
                        data.table = F, header = T, sep = ",")
acidic <- c("ASP", "GLU")
basic <- c("ARG", "LYS", "HIS")

eq_acid <- function(pKa, pH) {
  return(-1 / (1 + 10 ^ (pKa - pH)))
}

eq_base <- function(pKa, pH) {
  return(1 / (1 + 10 ^ (pH - pKa)))
}

# Calculate net charge
df <- df %>%
  mutate(aa = sub("_.*$", "", Residue)) %>%
  filter(aa %in% c(acidic, basic)) %>%
  mutate(charge = ifelse(aa %in% acidic,
                         eq_acid(pKa, 7.5),
                         eq_base(pKa, 7.5))) %>%
  group_by(PDB_ID) %>%
  summarise(total_charge = sum(charge))
write.csv(df, file="./database_charge_pH7.5.csv", quote = F, row.names = F)

# Map to Ensembl gene ID
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
annot <- data.table::fread("../annotation/HUMAN_9606_idmapping.dat", sep = "\t",
                           stringsAsFactors = F, data.table = F, header = F)
colnames(annot) <- c("uniprot", "database", "ID")
uniprot_pdb <- annot %>%
  filter(database == "PDB") %>%
  rename(PDB_ID = ID) %>%
  select(PDB_ID, uniprot)
df <- left_join(df, uniprot_pdb, by = "PDB_ID")
uniprot_ensembl <- annot %>%
  filter(database == "Ensembl") %>%
  rename(ensembl = ID) %>%
  select(ensembl, uniprot)
df <- inner_join(df, uniprot_ensembl, by = "uniprot")

# Expression data and cancer stage information
load("../expression_FPKM/BRCA.RData")
annot <- data.table::fread("../expression_FPKM/annotation/annot.tsv", data.table = F,
                           stringsAsFactors = F, header = T)

transform.data <- function(df, annot, project) {
  #' Separate different stages, and melt for ggplot.
  #'
  #' Require reshape2 to work (uses melt function).
  #'  
  #' @param df The data frame whose rownames are to be converted.
  #' @param annot Stage information.
  #' @param project The TCGA project name.
  #'
  #' @return A melted data.frame with different stages.
  library(reshape2)
  # df <- df[rowSums(df == 0) <= (ncol(df) / 2), ]
  df$symbol <- rownames(df)
  df <- melt(df, id.vars = "symbol")
  df$stage <- 0
  stage1 <- annot$barcode[grepl("(^i$)|(\\si[abc]?$)|(1)", annot$tumor_stage) &
                            annot$project == project]
  stage2 <- annot$barcode[grepl("(^ii$)|(\\si{2}[abc]?$)|(2)", annot$tumor_stage) &
                            annot$project == project]
  stage3 <- annot$barcode[grepl("(^iii$)|(\\si{3}[abc]?$)|(3)", annot$tumor_stage) &
                            annot$project == project]
  stage4 <- annot$barcode[grepl("(^iv$)|(\\siv[abc]?$)|(4)", annot$tumor_stage) &
                            annot$project == project]
  df$stage[df$variable %in% stage1] <- "i"
  df$stage[df$variable %in% stage2] <- "ii"
  df$stage[df$variable %in% stage3] <- "iii"
  df$stage[df$variable %in% stage4] <- "iv"
  df <- df[df$stage != 0, ]
  return(df)
}

process_exp <- function(exp, charge, annot, sample) {
  require(reshape2)
  exp <- exp[sub("\\..*$", "", rownames(exp)) %in% charge$ensembl, ]
  if (sample == "normal") {
    # exp <- exp[rowSums(exp == 0) <= (ncol(exp) / 2), ]
    exp$symbol <- rownames(exp)
    exp <- melt(exp, id.vars = "symbol")
    exp$stage <- "control"
  }
  else {
    exp <- transform.data(exp, annot, "BRCA")
  }
  exp$symbol <- sub("\\..*$", "", exp$symbol)
  return(exp)
}

dfn <- process_exp(datan, df, annot, "normal")
dft <- process_exp(datat, df, annot, "tumor")
exp <- rbind.data.frame(dfn, dft)
exp <- exp %>%
  group_by(symbol, stage) %>%
  summarise(expression = mean(value)) %>%
  rename(ensembl = symbol)
df <- df %>%
  select(ensembl, total_charge) %>%
  group_by(ensembl) %>%
  summarise(net_charge = mean(total_charge))
df <- inner_join(df, exp, by = "ensembl")
df <- df %>%
  mutate(value = log2(expression + 0.1) * net_charge)

# Net charge distribution
x <- select(df, ensembl, net_charge)
x <- x[!duplicated(x$ensembl), ]
x$ensembl <- factor(x$ensembl, levels=x$ensembl[order(x$net_charge)])
chrg_dist <- ggplot(x, aes(ensembl, net_charge)) +
             geom_col() +
             theme(axis.text.x = element_blank()) +
             xlab("Gene") +
             ylab("Net Charge") +
             ggtitle("Protein net charge distribution")
ggsave(chrg_dist, filename = "./charge_distribution_pH7.5.tiff", device = "tiff",
       width = 16, height= 9, units = "in", dpi = 200)

# Expression distribution
x <- x[order(x$net_charge), ]
df$ensembl <- factor(df$ensembl, levels=x$ensembl)
exp_dist <- ggplot(df, aes(ensembl, expression, color=stage)) +
            geom_col() +
            geom_vline(xintercept = tail(which(abs(x$net_charge) <= 0.0001),1),
                       linetype = "dotted") +
            facet_wrap(~stage) +
            theme(axis.text.x = element_blank()) +
            xlab("Gene") +
            ylab("log2(FPKM + 0.1)") +
            ggtitle("Gene Expression Distribution")
ggsave(exp_dist, filename = "./expression_distribution.tiff", device = "tiff",
       width = 16, height= 9, units = "in", dpi = 200)
# Protein net charge * expression
exp_chrg <- ggplot(df, aes(value, color=stage)) +
            geom_density() +
            xlab("Net Charge * Expression")+
            ggtitle("Charge * Expression Distribution")
ggsave(exp_chrg, filename = "./exp_charge_pH7.5.tiff", device = "tiff",
       width = 12, height= 9, units = "in", dpi = 200)
