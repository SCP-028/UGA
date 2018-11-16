# If we know the protein content of a given cell, then protein mass can be
# calculated (knowing the # and molecular weight of AAs). Assuming everything
# else is water, we then know the number of H2O molecules.
# According to http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/,
# there're in total 3e6 * 3e4 = 9e10 Da proteins, ~
rm(list=ls())

setwd("C:/Users/yz73026/Desktop/protein_pka/KAUST/result")
require(tidyverse)
require(glue)
require(ggpubr)

PROTEIN_PERCENT <- 0.165
WATER_MOLECULAR_WEIGHT <- 18.01528
WATER_ION_PRODUCT <- 1.55e-7

df <- data.table::fread("./pka_dataframe.csv",
                        stringsAsFactors = F, data.table = F, header = T, sep = ",")
acidic <- c("ASP", "GLU")
basic <- c("ARG", "LYS", "HIS")

eq_acid <- function(pKa) {
  return(sqrt(10 ^ -pKa))
}

eq_base <- function(pKb) {
  return(-(10 ^ -pKb))
}

# calculate Ka or Kb value for each residue
df <- df %>%
  filter(Residue %in% c(acidic, basic)) %>%
  mutate(pH = ifelse(Residue %in% acidic,
                     eq_acid(pK),
                     eq_base(pK))) %>%
  group_by(Protein) %>%
  summarise(protein_pH = sum(pH))

# map to Ensembl gene ID
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
# we have to use uniprot ID as an intermediate
ensemblAnnot <- data.table::fread("C:/Users/yz73026/Desktop/annotation/HUMAN_9606_idmapping.dat",
                           sep = "\t", stringsAsFactors = F, data.table = F, header = F)
colnames(ensemblAnnot) <- c("uniprot", "database", "ID")
df <- ensemblAnnot %>%
  filter(database == "PDB") %>%
  rename(Protein = ID) %>%
  select(Protein, uniprot) %>%
  right_join(df, by="Protein")

uniprot_ensembl <- ensemblAnnot %>%
  filter(database == "Ensembl") %>%
  rename(ensembl = ID) %>%
  select(ensembl, uniprot)
df <- inner_join(df, uniprot_ensembl, by = "uniprot")
# save unique uniprot IDs
write.csv(unique(df$uniprot), file = "unique_uniprot.csv", row.names = F, quote = F)
# ph: data.frame w/ columns: ensembl & overall_pH
ph <- df %>%
  select(Protein, ensembl, protein_pH) %>%
    group_by(ensembl) %>%
    summarise(overall_pH=mean(protein_pH))

# https://web.expasy.org/cgi-bin/compute_pi/pi_tool to retrieve protein molecular weight
mw <- read.table("protein_mass.tsv")
colnames(mw) <- c("pdb", "uniprot", "pI", "weight")
# mw: data.frame w/ columns uniprot & weight
mw <- mw %>%
  group_by(uniprot) %>%
  summarise(mw=mean(weight)) %>%
  left_join(uniprot_ensembl, by="uniprot")


estimate_ph <- function(exp, ph, mw) {
  #' Estimate environment pH given gene expression information.
  #'
  #' Require tidyverse to work.
  #'
  #' @param exp Data.frame with columns ensembl & expression.
  #' @param ph Data.frame with columns ensembl & overall_pH.
  #' @param mw Data.frame with columns ensembl & mw (molecular weight).
  #'
  #' @return A pH value (might be way out of range).
  rownames(exp) <- sub("\\..*$", "", rownames(exp))
  exp <- rowMeans(exp[rownames(exp) %in% ph$ensembl, ])
  exp <- data.frame(ensembl=names(exp), expression=exp)
  mw <- mw %>%
    inner_join(exp, by="ensembl") %>%
    mutate(product=expression * mw) %>%
    select(ensembl, product)
  water <- sum(mw$product) * (1 - PROTEIN_PERCENT) / (PROTEIN_PERCENT * WATER_MOLECULAR_WEIGHT)
  ph <- ph %>%
    inner_join(exp, by="ensembl") %>%
    mutate(product= expression * overall_pH) %>%
    select(ensembl, product)
  ph <- 7 - log10(1 + sum(ph$product) / (WATER_ION_PRODUCT * water))
  return (ph)
}


# load expression data and cancer stage information
filenames <- list.files("C:/Users/yz73026/Desktop/expression_FPKM/", pattern=".RData$", full.names = T)
projects <- sub(".*/([A-Z-]+)\\.RData", "\\1", filenames)
for (i in seq_along(projects)) {
  load(filenames[i])  # datan, datat, annot
  phn <- estimate_ph(exp=datan, ph=ph, mw=mw)
  pht <- estimate_ph(exp=datat, ph=ph, mw=mw)
  result <- rbind.data.frame(dfn, dft)
}
