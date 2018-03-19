#!/home/yi/miniconda3/bin/python
import os
import re
import glob
os.chdir('/home/yi/data/DNA_Methylation/LINE1/Data/matrix')
for filepath in glob.iglob('./**/*.txt', recursive=True):
    filename = re.search(r"[^\/][\w\-_]*\.txt$", filepath).group(0)
    with open('%s' % filepath, 'r') as file, \
            open(os.path.join('~/data/DNA_Methylation/LINE1/Data/annot/%s'
                 % filename), 'w') as annot, \
            open(os.path.join('~/data/DNA_Methylation/LINE1/Data/methy/%s'
                 % filename), 'w') as methy:
        for row in file:
            line = re.split('\t|\n', row)
            if re.match('!.*Sample_geo_accession', line[0]) is not None:
                annot.write(row)
            if re.match('!.*Sample_title', line[0]) is not None:
                annot.write(row)
            if re.match('!.*Sample_source_name', line[0]) is not None:
                annot.write(row)
            if re.match('\"ID_REF\"', line[0]) is not None:
                methy.write(row)
            if re.match(r'[^!]((cg\d*)|(\d*)|(A_17_P\d*))',
                        line[0]) is not None:
                methy.write(row)
