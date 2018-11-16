#!/home/yi/miniconda3/bin/python
import os
import glob
import csv
import re
os.chdir("/lustre1/yz73026/Bisulfite/GDCdata")
projects = "TCGA-BLCA"
""" All other projects
"TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC", "TCGA-LUSC",
"TCGA-STAD", "TCGA-COAD", "TCGA-READ", "TCGA-GBM"
"""
annot = list(csv.reader(open('geneName.tabular'), delimiter='\t'))
# 'Gene_Name \t Refseq \t Chro \t Start \t End'
print(annot[0:4][2])
# ['A1BG-AS1', 'NR_015380', 'chr19', '58863335', '58866549']
# Read all filepaths:
for filepath in glob.iglob('./%s/**/*.bed' % projects, recursive=True):
    filename = re.search(r"[^\/][\w\.]*bed$", filepath).group(0)
    print('Running: %s' % filename)
    with open('%s' % filepath, 'r') as file, \
        open(os.path.join('/media/yi/Data/Bisulfite/New/%s/%s'
             % (projects, filename)), 'w') as newfile:
        next(file)  # Jump header line
        for row in file:
            row = re.split('\t|\n', row)
            # ['chr1', '395779', '395780', '.', '0', '.', '0.00', '1', '']
            for i in range(len(annot)):
                if (row[0] == annot[i][2] and
                   int(annot[i][3]) <=
                   int(row[1]) <
                   int(annot[i][4])):
                    newfile.write('%r\t%r\t%r\t%r\tBody\t%r\t%r\n'
                                  % (row[0], row[1], row[2],
                                     annot[i][0], row[6], row[7]))
                elif (row[0] == annot[i][2] and
                      (int(annot[i][3]) - 2000) <=
                      int(row[1]) <
                      int(annot[i][3])):
                    newfile.write('%r\t%r\t%r\t%r\tPromoter\t%r\t%r\n'
                                  % (row[0], row[1], row[2],
                                     annot[i][0], row[6], row[7]))
        print('Finished %s\!' % filename)
