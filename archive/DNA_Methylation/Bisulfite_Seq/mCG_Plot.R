setwd("project/xynlab/from_oflow/yizhou/GDCdata")
projects <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC", "TCGA-LUSC", "TCGA-STAD", "TCGA-COAD", "TCGA-READ")## "TCGA-GBM"
library(plyr)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(8)
annot <- read.csv("geneName.tabular", header = T, sep = "\t")## Created with Galaxy: https://usegalaxy.org/
annot$Promoter <- annot$Start-2000
foreach (i=1:length(projects)) %dopar%## File path reading stuff
{
	files <- list.files(paste("./", projects[i], "/", sep = ""), pattern = "\\.(bed)$", recursive = T)
	barcode <- read.csv(paste("./", projects[i], "/barcode.txt", sep = ""), header = T, sep = "\t")
	barcode[] <- lapply(barcode, as.character)
	files <- files[order(sub(".*\\/", "", files))]
	files <- paste("./", projects[i], "/", files, sep = "")
	barcode <- barcode[order(barcode[ ,2]), ]
	normalb <- numeric()
	normalp <- numeric()
	cancerb <- numeric()
	cancerp <- numeric()
	for (j in 1:length(files))## Read one file
	{
		temp <- data.table::fread(files[j], header = F, data.table = F)
		temp <- temp[!(is.na(temp$V7) | temp$V7 == "."), ]
		colnames(temp) <- c("Chro","Start","End","Name","Score","Strand","Methy_Percent","numCTreads")
		for (x in 1:nrow(temp))
		{
			for (y in 1:nrow(annot))
			{
				if(temp$Chro[x] == annot$Chro[y] && findInterval(temp$Start[x], c(annot$Start[y], annot$End[y])) == 1L)## [CG] inside gene body
				{
					temp$Group[x] <- paste(temp$Group[x], "Body", sep = ";")
				}
				if(temp$Chro[x] == annot$Chro[y] && findInterval(temp$Start[x], c(annot$Promoter[y], annot$Start[y])) == 1L)## [CG] within upstream TSS2000
				{
					temp$Group[x] <- paste(temp$Group[x], "Promoter", sep = ";")
				}
			}
		}
		bodypct <- temp[grep("Body", temp$Group, perl = T), ]
		prompct <- temp[grep("Promoter", temp$Group, perl = T), ]
		if (grepl("TCGA-..-....-0\\w{2}",barcode[j,1]))## Cancer
		{
			percent <- sum(as.numeric(bodypct$Methy_Percent) * bodypct$numCTreads * 0.01) / sum(bodypct$numCTreads)
			cancerb <- c(cancerb,percent)
			percent <- sum(as.numeric(prompct$Methy_Percent) * prompct$numCTreads * 0.01) / sum(prompct$numCTreads)
			cancerp <- c(cancerp,percent)
		}
		else## Normal
		{
			percent <- sum(as.numeric(bodypct$Methy_Percent) * bodypct$numCTreads * 0.01) / sum(bodypct$numCTreads)
			normalb <- c(normalb,percent)
			percent <- sum(as.numeric(prompct$Methy_Percent) * prompct$numCTreads * 0.01) / sum(prompct$numCTreads)
			normalp <- c(normalp,percent)
		}
	}
	normalp <- data.frame(Type = "Normal", Group = "Promoter", Value = as.numeric(normalp))
	normalb <- data.frame(Type = "Normal", Group = "Body", Value = as.numeric(normalb))
	cancerp <- data.frame(Type = "Cancer", Group = "Promoter", Value = as.numeric(cancerp))
	cancerb <- data.frame(Type = "Cancer", Group = "Body", Value = as.numeric(cancerb))

	plot.data <- rbind(normalp,normalb,cancerp,cancerb)
	ggplot(plot.data,aes(x = Type,y = Value,fill = Group))+
	geom_boxplot()+
	labs(title = projects[i])
	ggsave(paste(projects[i],".png"),path = "./plot/")
	write.table(plot.data,paste("./plot/",projects[i],sep = ""))
}