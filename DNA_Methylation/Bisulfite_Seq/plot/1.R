setwd("/lustre1/yz73026/Bisulfite/GDCdata/New")
# project <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC", "TCGA-LUSC", "TCGA-STAD", "TCGA-COAD", "TCGA-READ")  # TCGA-GBM
project <- "TCGA-BLCA"
library(plyr)
library(ggplot2)
files <- list.files(paste("./", project, "/", sep = ""), pattern = "\\.(bed|BED)$")
barcode <- read.csv(paste("../barcode/", project, "/barcode.txt", sep = ""), header = T, sep = "\t")
barcode[] <- lapply(barcode, as.character)
files <- files[order(sub(".*\\/", "", files))]
files <- paste("./", project, "/", files, sep = "")
barcode <- barcode[order(barcode[ ,2]), ]
normalb <- numeric()
normalp <- numeric()
tumorb <- numeric()
tumorp <- numeric()
mCGnb <- numeric()
mCGnp <- numeric()
mCGtb <- numeric()
mCGtp <- numeric()
for (i in 1:length(files))
{
	temp <- data.table::fread(files[i], header = F, data.table = F)
	colnames(temp) <- c("Chromosome", "Start", "End", "Gene_Symbol", "Gene_Group", "Methy_Percent", "numCTreads")
	temp <- temp[!(temp$Methy_Percent == "'.'"), ]
	temp$Methy_Percent <- gsub("'", "", temp$Methy_Percent)
	temp$numCTreads <- gsub("'", "", temp$numCTreads)
	temp$mCG <- round(as.numeric(temp$Methy_Percent) * as.numeric(temp$numCTreads) * 0.01)
	tempb <- temp[temp$Gene_Group == "Body", ]
	tempp <- temp[temp$Gene_Group == "Promoter", ]
	rm(temp)  # Save memory
	percentb <- mean(as.numeric(tempb$Methy_Percent))
	percentp <- mean(as.numeric(tempp$Methy_Percent))
	mCGb <- mean(as.numeric(tempb$mCG))
	mCGp <- mean(as.numeric(tempp$mCG))
	if (grepl("TCGA-..-....-0\\w{2}", barcode[i,1]))
	{
		tumorb <- c(tumorb, percentb)
		tumorp <- c(tumorp, percentp)
		mCGtb <- c(mCGtb, mCGb)
		mCGtp <- c(mCGtp, mCGp)
	}
	if (grepl("TCGA-..-....-[12]\\w{2}", barcode[i,1]))
	{
		normalb <- c(normalb, percentb)
		normalp <- c(normalp, percentp)
		mCGnb <- c(mCGnb, mCGb)
		mCGnp <- c(mCGnp, mCGp)
	}
}
rm(tempb, tempp)
normalb <- data.frame(Type = "Normal", Group = "Body", Methy_Percent = as.numeric(normalb), mCG = as.numeric(mCGnb))
normalp <- data.frame(Type = "Normal", Group = "Promoter", Methy_Percent = as.numeric(normalp), mCG = as.numeric(mCGnp))
tumorb <- data.frame(Type = "Cancer", Group = "Body", Methy_Percent = as.numeric(tumorb), mCG = as.numeric(mCGtb))
tumorp <- data.frame(Type = "Cancer", Group = "Promoter", Methy_Percent = as.numeric(tumorp), mCG = as.numeric(mCGtp))
plot.data <- rbind(normalp,normalb,tumorp,tumorb)

ggplot(plot.data, aes(x = Type, y = Methy_Percent, fill = Group))+
geom_boxplot()+
labs(title = paste(project, "Methy_Percent", sep = "_"))
ggsave(paste(project, "_Methy_Percent", ".png", sep = ""), path = "./plot")

ggplot(plot.data, aes(x = Type, y = mCG, fill = Group))+
geom_boxplot()+
labs(title = paste(project, "mCG", sep = "_"))
ggsave(paste(project, "_mCG", ".png", sep = ""), path = "./plot")

write.table(plot.data, paste("./plot", project, sep = ""))
