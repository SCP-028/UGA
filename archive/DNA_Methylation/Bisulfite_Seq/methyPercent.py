#!python3
import os
import glob
import re
import csv
os.chdir("/lustre1/yz73026/Bisulfite/GDCdata")
project = "TCGA-LUSC"
barcode = list(csv.reader(open("./barcode/%s/barcode.txt" % project),
                          delimiter="\t"))
for filepath in glob.iglob("./%s/*.bed" % project,
                           recursive=True):
    filename = re.search(r"[^\/][\w\.]*bed$", filepath).group(0)
