""" Used after cor_coexpression.R """
import pandas as pd
import glob
import re
import os

filenames = glob.glob("./data/aa_paper/*.tsv")
projects = set([re.sub(r"^.*/(\w+)_\w+_gsea.*$", r"\1", f) for f in filenames])
pathways = set([re.sub(r"^.*/\w+_(\w+)_gsea.*$", r"\1", f) for f in filenames])

os.chdir("./aa_paper/")

for pw in pathways:
    writer = pd.ExcelWriter(f'{pw}.xlsx', engine='xlsxwriter')
    for proj in projects:
        f = f"../data/aa_paper/{proj}_{pw}_gsea.tsv"
        df = pd.read_table(f)
        df.to_excel(writer, sheet_name=proj, index=False)
    writer.save()
    writer.close()
