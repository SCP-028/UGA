# Visualization of data from Prof. Ying Li from Jilin University
rm(list=ls())
library(ggpubr)
library(glue)
setwd("C:/Users/yz73026/Desktop")
filenames <- list.files("./aa_dist", pattern="Rdata$", full.names=T)
projects <- toupper(sub(".*_(\\w+)_tcga.*$", "\\1", filenames))
for (i in seq_along(filenames)) {
    load(filenames[i])
    df <- data.frame(aa=names(ratio_tumor_normal), value=ratio_tumor_normal)
    # df$value[df$value < 1] <- -1 / df$value[df$value < 1]
    df$value <- log2(df$value + 0.1)
    df$group <- ifelse(df$value > 0, "tumor", "normal")
    df$aa <- factor(df$aa, levels=df$aa[order(df$value, decreasing=T)])
    p <- ggplot(df, aes(x=aa, y=value))+
         geom_bar(stat="identity", width=0.5)+
         theme_minimal()+
         ggtitle(projects[i])+
         xlab("")+
         ylab("log2(Tumor / Normal)")+
         theme(
             axis.ticks=element_blank(),
             axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.3, size=12),
             plot.margin=grid::unit(c(0,0,0,0), "mm")
             )
    p <- ggbarplot(
        df, x="aa", y="value", fill="group", color="white", palette="jco",
        x.text.angle=90, legend="none", ggtheme=theme_light(),
        title=projects[i], xlab=F, ylab="log2(Tumor / Normal)"
        )
    ggsave(filename=glue("{projects[i]}.tiff"),
         plot=p, device="tiff", width=6, height=4, units="in", dpi=300)
}

# sessionInfo()
# R version 3.4.4 (2018-03-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)

# Matrix products: default

# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base

# other attached packages:
# [1] glue_1.2.0    ggpubr_0.1.6  magrittr_1.5  ggplot2_2.2.1

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.16     assertthat_0.2.0 dplyr_0.7.4      grid_3.4.4       plyr_1.8.4       R6_2.2.2         gtable_0.2.0     scales_0.5.0     pillar_1.2.1     rlang_0.2.0      lazyeval_0.2.1   bindrcpp_0.2.2
# [13] labeling_0.3     ggsci_2.8        tools_3.4.4      purrr_0.2.4      munsell_0.4.3    compiler_3.4.4   pkgconfig_2.0.1  colorspace_1.3-2 bindr_0.1.1      tibble_1.4.2
