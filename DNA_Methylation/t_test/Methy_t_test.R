####################################
##  t.test of datan_M and datat_M ##
####################################
setwd("/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData")
load("~/data/DNA_Methylation/Methy_ZhouYi/RData/COAD.TNMethylation450.RData")
datat_NA <- datat[complete.cases(datat),]
datan_NA <- datan[complete.cases(datan),]
datat_final <- datat_NA[rownames(datat_NA) %in% rownames(datan_NA), ]
datan_final <- datan_NA[rownames(datan_NA) %in% rownames(datat_final), ]
datat_M <- log2(datat_final/(1-datat_final))
datan_M <- log2(datan_final/(1-datan_final))
rm(datat,datan,datat_NA,datan_NA)
## 375,012 genes in data*_final, and 388,134 entries in datan_NA
## data*_M is the M-value of the samples

test_result <- matrix(nrow = nrow(datan_M),ncol = 4)# Create result matrix
colnames(test_result) <- c("Normal","Cancer","p_t","Methylation")
rownames(test_result) <- rownames(datan_M)

mean_n <- rowMeans(datan_final)
mean_t <- rowMeans(datat_final)# Calculate row means for later comparison

for(i in 1:375012)
{
	test_result[i,1] <- mean_n[i]
	test_result[i,2] <- mean_t[i]	
	Var <- var.test(datan_M[i,],datat_M[i,])$p.value
	if(Var <= 0.05)
	{
		test_result[i,3] <- t.test(datan_M[i,],datat_M[i,],var.equal = F, paired = F)$p.value
	}
	else
	{
		test_result[i,3] <- t.test(datan_M[i,],datat_M[i,],var.equal = T, paired = F)$p.value
	}
	if(test_result[i,3] <= 0.05)
	{
		if(mean_t[i] > mean_n[i])
		{
			test_result[i,4] <- "hyper"
		}
		else
		{
			test_result[i,4] <- "hypo"
		}
	}
	else
	{
		test_result[i,4] <- "-"
	}
	print(i)
}
write.csv(test_result,"/home/yi/data/DNA_Methylation/Methy_ZhouYi/t_test_result.txt")
