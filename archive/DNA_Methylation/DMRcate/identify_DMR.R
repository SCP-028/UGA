## Following https://www.bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf
## source("https://bioconductor.org/biocLite.R")
## biocLite("DMRcate", dependencies = T)
load("~/data/DNA_Methylation/Methy_ZhouYi/RData/COAD.TNMethylation450.RData")# load original data
datat_NA <- datat[complete.cases(datat),]
datan_NA <- datan[complete.cases(datan),]# remove rows that contain NA
datat_final <- datat_NA[rownames(datat_NA) %in% rownames(datan_NA), ]
datan_final <- datan_NA[rownames(datan_NA) %in% rownames(datat_final), ]# only keep those existing in both matrices
colnames(datan_final) <- paste(colnames(datan_final),"Normal",sep = "_")
colnames(datat_final) <- paste(colnames(datat_final),"Cancer",sep = "_")# add name-tag to matrices
data_beta <- cbind(datat_final,datan_final)# merge into 1 matrix

library(DMRcate)# lib
data_M <- log2(data_beta/(1-data_beta))# beta-value to M-value
data("dmrcatedata")# import potential SNP-confounded probes database
data_M.noSNPs <- rmSNPandCH(data_M,dist = 2,mafcut = 0.05)# remove SNPs, filter out probes 2 nucleotides or closer toa SNP that have a minor allele frequency greater than 0.05, and the approximately 30,000 cross-reactive probes, so as to reduce confounding.
patient <- factor(sub("_.*","",colnames(data_M)))# cases
type <- factor(sub(".*_","",colnames(data_M)))# Normal/Cancer
design <- model.matrix(~type)# set up design matrix
COAD <- subset(data_M.noSNPs, grepl("(^cg|^ch)",rownames(data_M.noSNPs)))# remove SNP probes in matrix
myannotation <- cpg.annotate("array", COAD,annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"), analysis.type="differential",design=design, coef=2)# uses ilmn12.hg19 by default, can change to hg38 and compare; coef: see page https://support.bioconductor.org/p/74759/
## Your contrast returned 173565 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
dmrcoutput <- dmrcate(myannotation,lambda = 1000, C = 2)# find most differentially methylated regions. lambda is set at default 1000bp
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")# convert DMR list to a GRanges object

result.ranges
## GRanges object with 24994 ranges and 6 metadata columns:
## 	seqnames	ranges 	strand |   no.cpgs	minfdr  Stouffer	maxbetafc	meanbetafc
## 	<Rle>		<IRanges>	<Rle> | <integer>	<numeric> <numeric>	<numeric>	<numeric>
## chr6:33129988-33153499	chr6	[33129988, 33153499]	* |       152	4.442258e-227	0	0.3550392	0.1261248
## chr11:31817810-31841980	chr11	[31817810, 31841980]	* |       103	9.130968e-262	0	-0.4676111	-0.1933964
## chr2:63273436-63287288	chr2	[63273436, 63287288]	* |        90	0.000000e+00	0	-0.6186774	-0.1675779
## chr11:2158555-2165961	chr11	[ 2158555,  2165961]	* |        73	0.000000e+00	0	-0.4038374	-0.1674818
## chr20:57424521-57431303	chr20	[57424521, 57431303]	* |        72	4.674327e-279	0	-0.2508092	-0.1068670
##                        ...      ...
## chr7:5125074-5125232	chr7	[  5125074,   5125232]	* |         2	0.9783979	0.9898077	-0.0076106500	-0.0047949761
## chr14:21673491-21673796	chr14	[ 21673491,  21673796]	* |         2	0.9957989	0.9926927	0.0045255412	0.0028298160
## chrX:153084485-153085101	chrX	[153084485, 153085101]	* |         5	0.5577587	0.9933143	-0.0806472051	-0.0232543172
## chr7:29725145-29725196	chr7	[ 29725145,  29725196]	* |         3	0.9281623	0.9945085	-0.0003941554	-0.0002434241
## chrX:122865777-122867073	chrX	[122865777, 122867073]	* |         6	0.8606529	0.9958869	0.0264287047	0.0106980007
## overlapping.promoters
## <character>
## chr6:33129988-33153499	COL11A2-004
## chr11:31817810-31841980	PAX6-002, PAX6-011, PAX6-004, PAX6-005, PAX6-001, PAX6-003, PAX6-009, PAX6-008, PAX6-201, PAX6-012, PAX6-013, PAX6-016, PAX6-010, PAX6-006, PAX6-014, PAX6-015, PAX6-007, RCN1-011, PAX6-021, PAX6-020, PAX6-022, RCN1-013, PAX6-028, PAX6-025, PAX6-027, PAX6-023, RCN1-001, PAX6-029, PAX6-019, PAX6-024, PAX6-026
## chr2:63273436-63287288	OTX1-001, OTX1-002, OTX1-005, AC009501.4-004, AC009501.4-003, AC009501.4-001, AC009501.4-002, EHBP1-019, OTX1-006, OTX1-004, OTX1-003
## chr11:2158555-2165961	IGF2-AS-001, IGF2-AS-003, IGF2-009, IGF2-004, IGF2-003, IGF2-002, IGF2-001, IGF2-005, IGF2-007, IGF2-AS-002
## chr20:57424521-57431303	GNAS-050, GNAS-037, GNAS-058, GNAS-001, GNAS-009, GNAS-049, GNAS-AS1-001, GNAS-036, GNAS-059, GNAS-012, GNAS-AS1-004
## chr7:5125074-5125232	<NA>
## chr14:21673491-21673796	LINC00641-002, LINC00641-001
## chrX:153084485-153085101	<NA>
## chr7:29725145-29725196	ZNRF2P2-003, ZNRF2P2-002, ZNRF2P2-001, DPY19L2P3-002, DPY19L2P3-003
## chrX:122865777-122867073	THOC2-201, THOC2-001, THOC2-017
## -------
##   seqinfo: 23 sequences from an unspecified genome; no seqlengths


########
# Plot #
########
#groups <- c(Cancer = "magenta", Normal = "forestgreen")
#cols <- groups[as.character(type)]# colors denoting phenotypes
#samps <- c(1:290,291:328)
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=data_beta, phen.col=cols, genome="hg19", samps=samps)# plot significant DMR

df <- data.frame(seqnames=seqnames(results.ranges),
					starts=start(results.ranges)-1,
					ends=end(results.ranges),
					width=ranges(results.ranges)@width,
					scores=(elementMetadata(results.ranges)$no.cpgs)
				)# chromosome, start bp, end bp, length, score. Score selected as the number of CpGs in the sequence. This is done following "Detection of differentially methylated regions in whole genome bisulfite sequencing data using local Getis-Ord statistics, Wen et al., 2016"
write.table(df, file="DMR.bed", quote=F, sep="\t", row.names=F, col.names=F)# save the BED file to query on GREAT. Note that the "Significant By Region-based Binomial" view should be used.
