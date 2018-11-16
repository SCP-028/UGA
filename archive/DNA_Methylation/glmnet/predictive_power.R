setwd("~/data/DNA_Methylation/Methy_array/RData/Glmnet/")
library(glmnet)
## Load gene name, annot file, and raw data
load("./Raw_Data/normal.RData")
load("./Raw_Data/cancer.RData")
namesn <- rownames(datan)
namest <- rownames(datat)
## Load normal cv.glmnet data
load("./GLM_normal.RData")
cpgn <- GLMdata
## Load cancer cv.glmnet data
load("./GLM_cancer.RData")
cpgt <- GLMdata
rm(GLMdata)
## R^2 storage data.frame
r2df <- data.frame(namest,rep(0,length(namest)))
colnames(r2df) <- c("Gene_Name","R_Square")
library(plyr)
## Group count data.frame
countGroup <- data.frame("1stExon","5'UTR","Body","TSS1500","TSS200","3'UTR")
countGroup <- as.data.frame(cbind(t(countGroup),rep(0,6)))
rownames(countGroup) <- 1:6
colnames(countGroup) <- c("Gene_Group","freq")
countGroup[ ,2] <- c(0,0,0,0,0,0)
## Loop over GLMs
for(i in 1:length(cpgt))
{
	if(class(cpgt[[i]]) == "matrix")## All these with more than 1 CpG islands
	{
		coef <- cpgt[[i]]
		coef <- as.data.frame(coef[(coef[ ,1] != 0), ,drop = F]) #  Coefficient non-0 CpG islands
		if(nrow(coef) > 1)## Actually have CpGs rather than just intercept
		{
			annoti <- annot[annot$IlmnID %in% rownames(coef), ]
			annoti <- annoti[annoti$Gene_Name == namest[i], ]## Annotation rows in GLM and correspond to the gene
			y <- datat[rownames(datat) == namest[i], ,drop = F]## Gene expression data
			y <- as.data.frame(t(y[!(is.na(y[,1]) | y[,1]==""), ,drop = F]))## Somehow creates lots of NA rows, remove them
			colnames(y) <- "y"## Easier for later use
			x <- as.data.frame(t(methyt[rownames(methyt) %in% annoti$IlmnID, ]))## Model with CpG islands that have !0 values
			df <- as.data.frame(cbind(y,x))## Formula data.frame
			fit <- lm(y ~ .,data = df)## Fit linear model
			r2df[i,2] <- summary(fit)$r.squared## Save R^2
			if(r2df[i,2] >= 0.7)## Select R^2 >= 0.7
			{
				annoti <- count(annoti,"Gene_Group")## Count occurence of different gene groups
				countGroup <- rbind(countGroup,annoti)## Save in one data.frame
				countGroup <- aggregate(freq~Gene_Group,data=countGroup,FUN=sum)
			}
		}		
	}
}
#write.table(r2df,"./GLM/R2_cancer.txt")
write.table(countGroup,"./GLM/0.70_Group_cancer.txt")