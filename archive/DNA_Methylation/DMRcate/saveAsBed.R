df <- data.frame(seqnames=seqnames(results.ranges),
                starts=start(results.ranges)-1,
                ends=end(results.ranges),
                names=c(rep(".", length(results.ranges))),
                scores=(elementMetadata(results.ranges)$meanbetafc),
                strands=strand(results.ranges))
write.table(df, file="DMR.bed", quote=F, sep="\t", row.names=F, col.names=F)