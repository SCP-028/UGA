#!/home/yi/miniconda3/bin/python
import os
import glob
import csv
import re
# os.chdir("/lustre1/yz73026/Bisulfite/GDCdata")
os.chdir("/media/yi/Data/Bisulfite/")

project = "TCGA-COAD"
""" All the projects:
"TCGA-BLCA","TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC",
"TCGA-LUSC", "TCGA-STAD", "TCGA-COAD", "TCGA-READ", "TCGA-GBM"
"""

annot = list(csv.reader(open('geneName.tabular'), delimiter='\t'))
chros = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
         "chr17", "chr17_gl000205_random""chr18", "chr19",
         "chr19_gl000209_random", "chr1_gl000192_random", "chr2", "chr20",
         "chr21", "chr22", "chr3", "chr4", "chr4_gl000193_random",
         "chr4_gl000194_random", "chr5", "chr6", "chr7",
         "chr7_gl000195_random", "chr8", "chr9", "chrM", "chrUn_gl000211",
         "chrUn_gl000212", "chrUn_gl000213", "chrUn_gl000215",
         "chrUn_gl000218", "chrUn_gl000219", "chrUn_gl000220",
         "chrUn_gl000222", "chrUn_gl000223", "chrUn_gl000224",
         "chrUn_gl000228", "chrUn_gl000241", "chrX", "chrY"]
chrodict = {}
for chro in chros:
    # Sort annot by different chromosomes and into gene and promoter group
    chrodict['%s' % chro] = []
    chrodict['promoter_%s' % chro] = []
    for line in annot:
        if line[2] == chro:
            newline = []
            chrodict['%s' % chro].append(line)
            newline = [line[0], line[1], line[2],
                       int(line[3]) - 2000, int(line[3])]
            chrodict['promoter_%s' % chro].append(newline)
# print(chrodict['promoter_chr1'][0])
# print(chrodict['chr1'][0])

for filepath in glob.iglob('./%s/**/*.bed' % project, recursive=True):
    filename = re.search(r"[^\/][\w\.]*bed$", filepath).group(0)
    # print('Running: %s' % filename)
    """ Used in Sapelo.
    with open('%s' % filepath, 'r') as file,
        open(os.path.join('/lustre1/yz73026/Bisulfite/GDCdata/Split/%s/%s'
             % (project, filename)), 'w') as newfile:
    """
    with open('%s' % filepath, 'r') as file, \
        open(os.path.join('/media/yi/Data/Bisulfite/Split/%s/%s'
             % (project, filename)), 'w') as newfile:
        next(file)  # Jump header line
        for row in file:
            row = re.split('\t|\n', row)
            for chro in chros:
                if row[0] == chro:
                    for i in range(len(chrodict['%s' % chro])):
                        if int(chrodict['%s' % chro][i][3]) \
                                <= int(row[1]) \
                                < int(chrodict['%s' % chro][i][4]):
                            newfile.write('%r\t%r\t%r\t%r\tBody\t%r\t%r\n'
                                          % (row[0], row[1], row[2],
                                             chrodict['%s' % chro][i][0],
                                             row[6], row[7]))
                        if chrodict['promoter_%s' % chro][i][3] \
                                <= int(row[1]) \
                                < chrodict['promoter_%s' % chro][i][4]:
                            newfile.write('%r\t%r\t%r\t%r\tPromoter\t%r\t%r\n'
                                          % (row[0], row[1], row[2],
                                             chrodict['promoter_%s'
                                                      % chro][i][0],
                                             row[6], row[7]))
                    # print('Finished %s\!' % filename)
