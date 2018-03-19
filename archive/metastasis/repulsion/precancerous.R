try(setwd("C:/Users/jzhou/Desktop/GPL570"), silent=T)
try(setwd("/home/yizhou/data/GPL570"), silent=T)

annot <- data.table::fread("GPL570_annot.csv", stringsAsFactors = F,
                           data.table = F,
                           select = c("Probe Set ID", "Gene Symbol",
                                      "Entrez Gene", "Gene Title",
                                      "Gene Ontology Biological Process"))
colnames(annot) <- c("probe", "symbol", "entrez", "description", "GOBP")
