setwd("/lustre1/yz73026/Glmnet/alt/")
load("./normal.RData")
library(glmnet)
library(foreach)
library(doMC)
registerDoMC(24)
GLMdata <- foreach(i =1:nrow(datan)) %dopar%
{
	y <- as.matrix(datan[i, ])
	colnames(y) <- rownames(datan)[i]
	z <- as.character(annot[(annot[ ,2] == colnames(y)),1])
	x <- as.matrix(methyn[rownames(methyn) %in% z, ])
	if(nrow(x) == 1 | ncol(x) == 1) cor.test(x,y)
	else if(nrow(x) >= 2)
	{
		x <- t(x)
		class(x) <- "numeric"
		cvfit <- cv.glmnet(x,y)
		as.matrix(coef(cvfit, s = "lambda.min"))
	}
	else 0
}
rm(x,y,methyn,datan,annot,cvfit)
save.image("./GLM_normal.RData")
