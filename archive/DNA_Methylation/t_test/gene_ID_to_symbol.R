#geneID 2 symbol
setwd("/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/piano")
load("~/data/DNA_Methylation/Methy_ZhouYi/RData/piano/gene_expression.RData")
gene_info <- read.csv("./Homo_sapiens.gene_info",sep = "\t",header = T)
gene_info <- gene_info[ ,2:3]
output <- character(nrow(compare))
entrez <- as.character(compare[ ,2])
symbol <- as.character(gene_info[ ,2])
geneID <- character(0)
condition <- grepl(";",entrez)
pb = txtProgressBar(min = 0, max = nrow(compare), initial = 0) 
for (i in 1:nrow(compare))
{
	if(condition[i])
	{
		geneID <- unlist(strsplit(entrez[i],";"))
		for (j in 1:length(geneID))
		{
			k <- match(geneID[j], gene_info[,1])
			geneID[j] <- symbol[k]
		}
		output[i] <- paste(geneID, collapse = ";")
		geneID <- character(0)
	}
	else
	{
		k <- match(compare[i,2], gene_info[,1])
		output[i] <- symbol[k]
	}
	setTxtProgressBar(pb,i)
}
compare <- cbind(compare,output)
colnames(compare) <- c("Refseq_mRNA","Entrez","Symbol")
write.csv(compare,"./Entrez_symbol.txt")