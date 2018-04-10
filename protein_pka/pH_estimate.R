# If we know the protein content of a given cell, then protein mass can be
# calculated (knowing the # and molecular weight of AAs). Assuming everything
# else is water, we then know the number of H2O molecules.
# According to http://book.bionumbers.org/what-is-the-macromolecular-composition-of-the-cell/,
# and 55% percent of the cell dry mass are proteins, water makes up 70% of the total cell weight,
# which means 16.5% of the total weight of cells are proteins.
rm(list=ls())

setwd("C:/Users/yz73026/Desktop/protein_pka/KAUST/result")
require(dplyr)
require(glue)
require(ggpubr)

df <- data.table::fread("./pka_dataframe.csv",
                        stringsAsFactors = F, data.table = F, header = T, sep = ",")
acidic <- c("ASP", "GLU")
basic <- c("ARG", "LYS", "HIS")

eq_acid <- function(pKa) {
  return(sqrt(10 ^ pKa))
}

eq_base <- function(pKb) {
  return(-(10 ^ pKb))
}

df <- df %>%
  filter(Residue %in% c(acidic, basic)) %>%
  mutate(pH = ifelse(Residue %in% acidic,
                     eq_acid(pK),
                     eq_base(pK))) %>%
  group_by(Protein) %>%
  summarise(protein_pH = sum(pH))

# Map to Ensembl gene ID
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
ensemblAnnot <- data.table::fread("C:/Users/yz73026/Desktop/annotation/HUMAN_9606_idmapping.dat",
                           sep = "\t", stringsAsFactors = F, data.table = F, header = F)
colnames(ensemblAnnot) <- c("uniprot", "database", "ID")
uniprot_pdb <- ensemblAnnot %>%
  filter(database == "PDB") %>%
  rename(Protein = ID) %>%
  select(Protein, uniprot)
df <- left_join(df, uniprot_pdb, by = "Protein")
uniprot_ensembl <- ensemblAnnot %>%
  filter(database == "Ensembl") %>%
  rename(ensembl = ID) %>%
  select(ensembl, uniprot)
df <- inner_join(df, uniprot_ensembl, by = "uniprot")
df <- df %>%
  select(ensembl, protein_pH) %>%
    group_by(ensembl) %>%
    summarise(overall_pH=sum(protein_pH))
# write.csv(df, file="./ensembl_charge_pH7.0.csv", quote=F, row.names=F)

# Expression data and cancer stage information
filenames <- list.files("C:/Users/yz73026/Desktop/expression_FPKM/", pattern=".RData$", full.names = T)
projects <- sub(".*/([A-Z-]+)\\.RData", "\\1", filenames)
for (i in seq_along(projects)) {
  load(filenames[i])  # datan, datat, annot
  rownames(datan) <- rownames(datat) <- sub("\\..*$", "", rownames(datan))
  datan <- rowMeans(datan[rownames(datan) %in% df$ensembl, ])
  datat <- rowMeans(datat[rownames(datat) %in% df$ensembl, ])
  datan <- data.frame(ensembl=names(datan), expression=datan)
  datat <- data.frame(ensembl=names(datat), expression=datat)
  dfn <- inner_join(datan, df, by="ensembl")
  dft <- inner_join(datat, df, by="ensembl")
  dfn <- dfn %>%
    mutate(
      charge_expression=ifelse(Net_charge >= 0,
                               log2(overall_pH * expression + 0.1),
                               -log2(-overall_pH * expression + 0.1)
                               ),
      sample="control",
      project=projects[i]
    )
  dft <- dft %>%
    mutate(
      charge_expression=ifelse(Net_charge >= 0,
                               log2(Net_charge * expression + 0.1),
                               -log2(-Net_charge * expression + 0.1)
                               ),
      sample="tumor",
      project=projects[i]
    )
  result <- rbind.data.frame(dfn, dft)
  p <- ggdensity(result, x = "charge_expression", add="mean", fill="sample", color="sample")
  p <- ggpar(p, ggtheme=theme_minimal(), palette = c("#00AFBB", "#E7B800"),
             title=projects[i], xlab="log2 (Net charge * expression + 0.1)", ylab="Density")
  ggsave(filename=glue("{projects[i]}.tiff"),
           plot=p, device="tiff", width=10, height=6, units="in", dpi=200)
}
