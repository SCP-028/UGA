import glob
import json
import subprocess

import pandas as pd

import xmltodict as xd


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


df = pd.read_pickle("/data/protein_pka/result/pI.pkl")

# Get PDB annotation
# downloading files using PDB's RESTful API
url = "https://www.rcsb.org/pdb/rest/describePDB?structureId="
pdb_error = []
pdb_urls = [url + ",".join(pdb) for pdb in chunks(df.protein.unique(), 1000)]
with open("/data/protein_pka/result/pdb_meta/url.list", "w") as f:
    f.write("\n".join(pdb_urls))
"""
r = subprocess.run(["aria2c", "--input-file", "/data/protein_pka/result/pdb_meta/url.list",
                              "-o", "/data/protein_pka/result/pdb_meta/describePDB"])
"""

# merging downloaded files into one json
filenames = glob.iglob("/data/protein_pka/result/pdb_meta/describe*")
ans = {}
for filename in filenames:
    with open(filename, "r") as f:
        pdb = xd.parse(f.read())
    for node in pdb['PDBdescription']['PDB']:
        ans[node["@structureId"]] = node

with open("/data/protein_pka/result/pdb_meta/pdb_description.json", "w") as f:
    f.write(json.dumps(ans))

PDB = pd.read_json("/data/protein_pka/result/pdb_meta/pdb_description.json").T

# Filter list, only keep human proteins
df = pd.read_table("/data/annotation/taxonomy/pdb_source.csv", header=None, sep=",")
df.columns = ["PDB", "taxonomy"]
df = df[~df["taxonomy"].isna()]
df = df[df["taxonomy"].str.contains("HOMO SAPIENS")]

PDB = PDB[PDB["@structureId"].isin(df["PDB"])]

# Convert to Ensembl Gene ID
annotation = pd.read_csv("/data/annotation/uniprot_id_map/HUMAN_9606_idmapping.csv")
uniprot_pdb = annotation[(annotation["ID"].isin(PDB.index)) & (annotation["ID_type"] == "PDB")]
uniprot_ensembl = annotation[annotation["ID_type"] == "Ensembl"]
annot = pd.merge(uniprot_pdb, uniprot_ensembl, on="UniProtKB_AC")[["ID_x", "ID_y"]]
annot.columns = ["PDB", "Ensembl"]

PDB = pd.merge(PDB, annot, left_on=PDB["@structureId"], right_on=annot["PDB"])
PDB = PDB.drop("key_0", axis=1)

# Sort GTEx data by cell type -> organ
ntissue = pd.read_table("/data/annotation/tissue_expression/human_protein_atlas_normal.tsv")
ttissue = pd.read_table("/data/annotation/tissue_expression/human_protein_atlas_tumor.tsv")
