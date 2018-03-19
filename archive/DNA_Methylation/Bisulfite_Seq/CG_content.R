try(setwd("/home/yi/data/DNA_Methylation/LINE1/Data/map"), silent = T)
try(setwd("/lustre1/yz73026/LINE1"), silent = T)
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-LUAD","TCGA-UCEC","TCGA-LUSC","TCGA-STAD","TCGA-COAD","TCGA-READ", "TCGA-GBM")  # GBM has no normal, so plotted separately.
barcode <- read.csv("barcode.txt", sep = "\t", header = T)
barcode <- as.data.frame(sapply(barcode, function(x) as.character(x)))
library(ggplot2)
library(dplyr)

eachPlot <- function(df, col)
{
#	Plot either mCG or mCG/CG
	if (col != "mCG" && col != "percentage")
	{
		stop(paste("No corresponding columns in df named: ", col, sep = ""))
	}
	if (col == "mCG")
	{
		normal <- as.numeric(df$mCG[df$type == "normal"])
		tumor <- as.numeric(df$mCG[df$type == "tumor"])
	}
	else if (col == "percentage")
	{
		normal <- as.numeric(df$percentage[df$type == "normal"])
		tumor <- as.numeric(df$percentage[df$type == "tumor"])
	}

	tumor <- data.frame(Group = "Tumor", Value = tumor)
	if (length(normal) != 0)
#	TCGA-GBM has no normal
	{
		normal <- data.frame(Group = "Normal", Value = normal)
		plot.data <- rbind(normal,tumor)
	}
	else
	{
		plot.data <- tumor
	}
	ggplot(plot.data, aes(x = Group, y = Value, fill = Group))+
	geom_boxplot()+
	labs(title = paste(projects[i], col, sep = "_"))
}

result <- matrix(nrow = 0, ncol = 7)
for(i in seq_along(projects))
{
#	Make file data frame, containing barcode, file name & path, project, type, mCG and percentage
	files <- list.files(paste("./",projects[i],"/",sep = ""), pattern = "\\.(bed|BED)$", recursive = T)
	files <- data.frame(filePath = files, fileName = sub(".*\\/","",files), project = projects[i], type = NA, mCG = NA, percentage = NA)
	files$filePath <- paste("./", projects[i], "/", files$filePath, sep = "")
	files <- inner_join(barcode, files, by = "fileName")
	files$type <- ifelse(grepl("TCGA-..-....-0\\w{2}", files$barcode), "tumor", "normal")

	for (j in seq_along(files$filePath))
	{
#		Read file and calculate average mCG and mCG/CG numbers
		temp <- data.table::fread(files$filePath[j], header = T, data.table = F)
		temp <- temp[!(is.na(temp$percentMeth) | temp$percentMeth == "."), ]
		files$mCG[j] <- sum(as.numeric(temp$percentMeth)* 0.01)
		files$percentage[j] <- mean(as.numeric(temp$percentMeth) * 0.01)
	}
	eachPlot(files, col = "mCG")
	ggsave(paste(projects[i], "_mCG", ".png", sep = ""),path = "./plot/")
	eachPlot(files, col = "percentage")
	ggsave(paste(projects[i], "_percentage", ".png", sep = ""),path = "./plot/")
#	write.table(files, paste("./plot/",projects[i],sep = ""))
	result <- rbind(result, as.matrix(files))
}

finalPlot <- function(df, col)
{
#	The only difference is that it has to separate different cancer types
	df <- as.data.frame(df)
	if (col != "mCG" && col != "percentage") 
	{
		stop(paste("No corresponding columns in df named: ", col, sep = ""))
	}
	if (col == "mCG")
	{
		normal <- as.numeric(df$mCG[df$type == "normal"])
		tumor <- as.numeric(df$mCG[df$type == "tumor"])
	}
	else if (col == "percentage")
	{
		normal <- as.numeric(df$percentage[df$type == "normal"])
		tumor <- as.numeric(df$percentage[df$type == "tumor"])
	}
	normal <- data.frame(Group = "Normal", Value = normal)
	tumor <- data.frame(Group = df$project[df$type == "tumor"], Value = tumor)
	plot.data <- rbind(normal,tumor)
	ggplot(plot.data, aes(x = Group, y = Value, fill = Group))+
	geom_boxplot()+
	labs(title = paste(projects[i], col, sep = "_"))
}

write.table(result, "./plot/resultTable.txt")
finalPlot(result, col = "mCG")
ggsave("result_mCG.png", path = "./plot/")
finalPlot(result, col = "percentage")
ggsave("result_percentage.png", path = "./plot/")
