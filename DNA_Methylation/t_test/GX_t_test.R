load("~/data/DNA_Methylation/Methy_ZhouYi/RData/GX_Methy.RData")
rm(compare,datan_hypo,datan_hyper,datat_hyper,datat_hypo,hyper,hypo,a)
test_result <- matrix(nrow = nrow(gexpn),ncol = 4)# Create result matrix
colnames(test_result) <- c("Normal","Cancer","p_t","Gene Expression")
rownames(test_result) <- rownames(gexpn)

mean_n <- rowMeans(gexpn)
mean_t <- rowMeans(gexpt)# Calculate row means for later comparison

for(i in 1:nrow(gexpn))
{
	test_result[i,1] <- mean_n[i]
	test_result[i,2] <- mean_t[i]	
	Var <- var.test(gexpn[i,],gexpt[i,])$p.value
	if(Var <= 0.05)
	{
		test_result[i,3] <- t.test(gexpn[i,],gexpt[i,],var.equal = F, paired = F)$p.value
	}
	else
	{
		test_result[i,3] <- t.test(gexpn[i,],gexpt[i,],var.equal = T, paired = F)$p.value
	}
	if(test_result[i,3] <= 0.05)
	{
		if(mean_t[i] > mean_n[i])
		{
			test_result[i,4] <- "up-regulated"
		}
		else
		{
			test_result[i,4] <- "down-regulated"
		}
	}
	else
	{
		test_result[i,4] <- "-"
	}
	print(i/14240*100)
}
write.csv(test_result,"/home/yi/data/DNA_Methylation/Methy_ZhouYi/GX_t_result.txt")