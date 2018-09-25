#!/home/yi/miniconda3/bin/python
import os
import glob
import csv
import re
from bisect import bisect_left
os.chdir("/lustre1/yz73026/Bisulfite/GDCdata")
project = "TCGA-COAD"
# "TCGA-BLCA","TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC",
# "TCGA-LUSC", "TCGA-STAD", "TCGA-READ", "TCGA-GBM"
annot = list(csv.reader(open('/lustre1/yz73026/LINE1/LINE1.txt'),
                        delimiter='\t'))
chros = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
         "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
         "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY",
         "chr1_gl000191_random", "chr1_gl000192_random",
         "chr4_gl000193_random", "chr4_gl000194_random",
         "chr7_gl000195_random", "chr8_gl000196_random",
         "chr8_gl000197_random", "chr9_gl000198_random",
         "chr9_gl000199_random", "chr9_gl000200_random",
         "chr9_gl000201_random", "chr11_gl000202_random",
         "chr17_gl000203_random", "chr17_gl000204_random",
         "chr17_gl000205_random", "chr17_gl000206_random",
         "chr18_gl000207_random", "chr19_gl000208_random",
         "chr19_gl000209_random", "chr21_gl000210_random",
         "chrUn_gl000211", "chrUn_gl000212", "chrUn_gl000213",
         "chrUn_gl000214", "chrUn_gl000215", "chrUn_gl000216",
         "chrUn_gl000217", "chrUn_gl000218", "chrUn_gl000219",
         "chrUn_gl000220", "chrUn_gl000221", "chrUn_gl000222",
         "chrUn_gl000223", "chrUn_gl000224", "chrUn_gl000225",
         "chrUn_gl000226", "chrUn_gl000227", "chrUn_gl000228",
         "chrUn_gl000229", "chrUn_gl000230", "chrUn_gl000231",
         "chrUn_gl000232", "chrUn_gl000233", "chrUn_gl000234",
         "chrUn_gl000235", "chrUn_gl000236", "chrUn_gl000237",
         "chrUn_gl000238", "chrUn_gl000239", "chrUn_gl000240",
         "chrUn_gl000241", "chrUn_gl000242", "chrUn_gl000243",
         "chrUn_gl000244", "chrUn_gl000245", "chrUn_gl000246",
         "chrUn_gl000247", "chrUn_gl000248", "chrUn_gl000249", "chrM"]
annotDict = {}
for chro in chros:
    """ ***** Break annotation file by chromosomes. ***** """
    annotDict['%s' % chro] = []
    for line in annot:
        if line[0] == chro:
            annotDict['%s' % chro].append(line)


def binary_search(a, x, lo=0, hi=None):
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)  # don't walk off the end


filename = 'ResultCount_MERGING_1_NIC1254A45.hg19_rCRSchrm.fa.realign.mdups.recal.cpg.filtered.sort.CG.6plus2.fixed.bed'
filepath = glob.glob('./%s/**/%s' % (project, filename), recursive=True)
with open('%s' % filepath[0], 'r') as file, \
    open(os.path.join('/lustre1/yz73026/Bisulfite/GDCdata/Split/%s/%s'
                      % (project, filename)), 'w') as newfile:
    next(file)
    newfile.write('chrom\tchromStart\tchromEnd\t'
                  'name\tpercentMeth\tnumCTreads\n')
    for row in file:
        row = re.split('\t|\n', row)
        for chro in chros:
            if row[0] == chro:
                for i in range(len(annotDict['%s' % chro])):
                    if (int(annotDict['%s' % chro][i][3]) <=
                            int(row[1]) <
                            int(annotDict['%s' % chro][i][4])):
                        newfile.write('%r\t%r\t%r\t%r\t%r\t%r\n'
                                      % (row[0], row[1], row[2],
                                         annotDict['%s' % chro][i][3],
                                         row[6], row[7]))
                        break
