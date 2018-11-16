t_result <- read.csv("/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/t_test/t_test_result.txt")
hyper <- t_result[t_result[ ,"Methylation"] == "hyper", ]
hypo <- t_result[t_result[ ,"Methylation"] == "hypo", ]
annot <- read.csv("/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/cpg2annotation450",sep = "\t")
hyper_annot <- annot[annot[,1] %in% hyper[,1], ]
hypo_annot <- annot[annot[,1] %in% hypo[,1], ]

hyper_order<- hyper[order(hyper[,1]), ]
hyper_annot_order<- hyper_annot[order(hyper_annot[,1]), ]
hyper_annot <- cbind(hyper_order[,1:4],hyper_annot_order[2:6])

hypo_order<- hypo[order(hypo[,1]), ]
hypo_annot_order<- hypo_annot[order(hypo_annot[,1]), ]
hypo_annot <- cbind(hypo_order[,1:4],hypo_annot_order[2:6])

write.csv(hyper_annot,"/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/hyper_annot.txt")
write.csv(hypo_annot,"/home/yi/data/DNA_Methylation/Methy_ZhouYi/RData/hypo_annot.txt")