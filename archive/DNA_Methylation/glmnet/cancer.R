setwd("/lustre1/yz73026/Glmnet/alt/")
load("./cancer.RData")
library(glmnet)
library(foreach)
library(doMC)
registerDoMC(24)
GLMdata <- foreach(i =1:nrow(datat)) %dopar%
{
	y <- as.matrix(datat[i, ])
	colnames(y) <- rownames(datat)[i]
	x <- as.matrix(methyt[rownames(methyt) %in% annot[annot[ ,2] == colnames(y),1], ])
	if(nrow(x) >= 2)
	{
		x <- t(x)
		class(x) <- "numeric"
		cvfit <- cv.glmnet(x,y)
		z <- coef(cvfit, s= "lambda.min")
		as.matrix(coef(cvfit, s = "lambda.min"))
	}
}
rm(x,y,methyt,datat,annot,cvfit)
save.image("./GLM_cancer.RData")
