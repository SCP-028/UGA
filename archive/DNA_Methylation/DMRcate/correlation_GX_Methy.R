title <- "hyper_GXD_Body_all"
setwd("~/data/DNA_Methylation/Methy_ZhouYi/RData/Methy_GX/meta")
## Read files
df <- read.csv("./methy_GX/Body/hyper_GXD_Body",header = T)
GX <- read.csv("./GX/GXall",header = T)
load("./methy/Methyall.RData")
## Remove duplicates
df$X <- NULL
df <- subset(df, !duplicated(subset(df, select = c(IlmnID, Entrez))))
## Get correct rownames of GX
a <- GX[ ,1]
a <- sub(".*\\|","",a)
rownames(GX) <- a
GX$X <- NULL
## Get correct colnames of methy
a <- colnames(methy)
a <- sub("(TCGA).(..).(\\d{4}).*","\\1_\\2_\\3",a)
colnames(methy) <- a
## Merge duplicated rows into one, only keeping the row with smallest Methy_p
GX <- GX[rownames(GX) %in% df$Entrez, ]
df_dup <- df[duplicated(df$Entrez) | duplicated(df$Entrez,fromLast = T), ]
df_nodu <- df[!df$Entrez %in% df_dup$Entrez, ]
df_dup <- df_dup[order(df_dup$Methy_P), ]
df_dup <- df_dup[!duplicated(df_dup$Entrez), ]
df <- rbind(df_nodu,df_dup)
rm(df_dup,df_nodu)
## IlmnID to Entrez ID
methy <- methy[rownames(methy) %in% df$IlmnID, ]
df <- df[order(df$IlmnID), ]
methy <- methy[order(rownames(methy)), ]
rownames(methy) <- df$Entrez
GX <- GX[order(rownames(GX)), ]
methy <- methy[order(rownames(methy)), ]
## Collect paired cases
GX <- GX[ ,colnames(GX) %in% colnames(methy)]
library(limma)
methy <- as.data.frame(avearrays(methy))
methy <- methy[ ,colnames(methy) %in% colnames(GX)]
methy <- t(methy)
GX <- t(GX)
GX <- GX[order(rownames(GX)), ]
methy <- methy[order(rownames(methy)), ]
GX <- GX[ ,order(colnames(GX))]
methy <- methy[ ,order(colnames(methy))]
## Run if first run failed:
#df_dup <- df[duplicated(df$IlmnID) | duplicated(df$IlmnID,fromLast = T), ]
#methy_dup <- methy[ ,colnames(methy) %in% df_dup$IlmnID]
#methy <- cbind(methy,methy_dup)
#methy <- methy[ ,colnames(methy) %in% df$IlmnID]
#methy <- methy[ ,order(colnames(methy))]
#a <- df$Entrez
#colnames(methy) <- a
#methy <- methy[ ,order(colnames(methy))]

## The other alternative
#colnames(methy) == df$IlmnID
#methy <- cbind(methy,methy[ ,"cg04722914"])
#methy <- methy[ ,order(colnames(methy))]
#colnames(methy)[1] <- "cg04722914"
#methy <- methy[ ,order(colnames(methy))]
#df <- df[order(df$IlmnID), ]
#colnames(methy) <- df$Entrez
#methy <- methy[ ,order(colnames(methy))]

## Correlate methy and GX
a <- data.frame(0,0)
colnames(a) <- c("P","Rho")
for(i in 1:ncol(GX))
{
	b <- cor.test(methy[ ,i],GX[ ,i],method = "spearman")
	a[i,1:2] <- c(b$p.value,b$estimate)
}
rownames(a) <- colnames(GX)
write.csv(a,paste("../result/correlation/",title))
## Plot: the larger bins is, the more specific
library(ggplot2)
ggplot(data = a, aes(a$Rho))+
geom_histogram(aes(y = ..density..), bins = 15)+
scale_x_continuous(breaks = seq(-1, 1, by = 0.1))+
geom_density(col = 2)+
labs(title = title)+
labs(x="Rho", y = "Density")
ggsave(paste(title,".png"),path = "../result/correlation/")
#qplot(a$Rho,geom = "histogram",main = title, xlab = "Rho", bins = 20)
