setwd("~/data/DNA_Methylation/Methy_ZhouYi/RData/Glmnet/")
load("./Raw_Data/normal.RData")
normal <- read.table("./GLM/Group_normal.txt")
library(plyr)
total <- count(annot,"Gene_Group")
df <- cbind(total,normal$freq_0.7)
colnames(df) <- c("Gene_Group","Total_Num","Predictive_Num")
df$Non_Predictive <- df$Total_Num - df$Predictive_Num
group <- as.character(df$Gene_Group)
rm(annot,datan,methyn,GLMdata,normal,total)
result <- matrix(rep(0,6),nrow = 6)
dimnames(result) <- list(group,"P_cai")
for(i in 1:length(group))
{
	x1 <- df[df$Gene_Group == group[i],3]
	x2 <- df[df$Gene_Group == group[i],4]
	x3 <- sum(df$Predictive_Num) - x1
	x4 <- sum(df$Non_Predictive) - x2
	temp <- matrix(c(x1,x2,x3,x4),nrow = 2,byrow = F)
	result[i,1] <- fisher.test(temp)$p.value
}
df <- cbind(df,result)
df$Ratio <- df$Predictive_Num / df$Non_Predictive
write.table(df,"./GLM/Group_Chi.txt")