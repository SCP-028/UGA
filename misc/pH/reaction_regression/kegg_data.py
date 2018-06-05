import os
import sys

import requests

ROOTPATH = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(ROOTPATH)
ss = requests.Session()

# download reactions file from KEGG
if not os.path.exists(f"{ROOTPATH}/kegg_reactions.tsv"):
    r = ss.get("http://rest.kegg.jp/list/reaction")
    if r.ok:
        with open(f"{ROOTPATH}/kegg_reactions.tsv", "w") as f:
            f.write(r.text)
df = pd.read_table(f"{ROOTPATH}/kegg_reactions.tsv", headers=False)
