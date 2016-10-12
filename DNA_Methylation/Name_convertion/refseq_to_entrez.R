#!/usr/bin/Rscript
setwd("/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/piano")
REF <- read.csv("./refseq.txt",sep = "\t",header = T)
ENT <- read.csv("./Entrez.txt",sep = "\t",header = T)
ENT <- ENT[,1:2]
REF <- REF[order(REF[ ,2]), ]
ENT <- ENT[order(ENT[ ,1]), ]
output <- character(nrow(REF))
refgene <- numeric(0)
condition <- grepl(";",REF[,2])
pb = txtProgressBar(min = 0, max = nrow(REF), initial = 0) 
for (i in 1:nrow(REF))
{
	if (condition[i])
	{
		refseq <- unlist(strsplit(as.character(REF[i,2]),";"))
		for (j in 1:length(refseq))
		{
			k <- match(refseq[j], ENT[,1])
			refgene[j] <- ENT[k,2]
		}
		output[i] <- paste(refgene, collapse = ";")
		refgene <- numeric(0)
	}
	else
	{
		k <- match(REF[i,2], ENT[,1])
		output[i] <- ENT[k,2]
	}
	setTxtProgressBar(pb,i)
}
REF <- cbind(REF,output)
colnames(REF) <- c("IlmnID","Refseq_mRNA","Entrez")
write.csv(REF,"./refseq_entrez.txt")