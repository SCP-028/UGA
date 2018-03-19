setwd("~/data/DNA_Methylation/Methy_ZhouYi/RData/")
# Methy data colnames -> GX colnames
load("./COAD.TN_noNA.RData")
methyn <- datan_final
methyt <- datat_final
rm(datan_final,datat_final)
colnames(methyn) <- sub("(TCGA)-(..)-(\\d{4}).*","\\1_\\2_\\3",colnames(methyn))
colnames(methyt) <- sub("(TCGA)-(..)-(\\d{4}).*","\\1_\\2_\\3",colnames(methyt))
write.table(methyn,"./methyn",sep = "\t")
write.table(methyt,"./methyt",sep = "\t")