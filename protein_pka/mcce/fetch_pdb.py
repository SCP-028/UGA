#!python3
import glob
import json
import subprocess

import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline

import xmltodict as xd

offline.init_notebook_mode()


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

# Sort by subcellular localization
protein_atlas = pd.read_table("/data/annotation/localization/protein_atlas_subcellular_location.tsv")
uniprot = pd.read_table("/data/annotation/localization/uniprot_locations_all.tsv")
reliability = ["Enhanced", "Supported", "Approved"]
ans = pd.DataFrame(columns=[
    "PDB", "Ensembl", "Gene name", "Localization", "Reliability",
    "@title", "@expMethod", "@resolution", "@last_modification_date"
])
for item in reliability:
    tmp = protein_atlas[protein_atlas["Reliability"] == item][["Gene", "Gene name", item]]
    tmp = tmp.rename(columns={item: "Localization"})
    tmp = tmp[~tmp["Gene"].isin(ans["Ensembl"])]
    tmp = pd.merge(PDB, tmp, left_on="Ensembl", right_on="Gene", how="inner").drop("Gene", axis=1)
    tmp["Reliability"] = item
    ans = pd.concat([ans, tmp], ignore_index=True)

ans = ans[[
    "PDB", "Ensembl", "Gene name", "Localization", "Reliability",
    "@title", "@expMethod", "@resolution", "@last_modification_date"
]]
ans = ans.dropna()
ans.to_pickle("/data/protein_pka/result/PDB_description.pkl")
localization = dict([[key, []] for key in set([x for item in ans["Localization"].unique() for x in item.split(";")])])
for i, row in ans.iterrows():
    for organelle in row["Localization"].split(";"):
        localization[organelle] += [[row["PDB"], row["Ensembl"], row["Gene name"]]]

pI = pd.read_pickle("/data/protein_pka/result/pI.pkl")

ans = pd.DataFrame()
organelle = "Nucleoplasm"
tmp = pI[pI["protein"].isin([x[0] for x in localization[organelle]])]
tmp["Localization"] = organelle
ans = pd.concat([ans, tmp], ignore_index=True)

for organelle in localization.keys():
    tmp = pI[pI["protein"].isin([x[0] for x in localization[organelle]])]
    tmp["Localization"] = organelle
    ans = pd.concat([ans, tmp], ignore_index=True)

organelles = ans["Localization"].unique()
data = [go.Histogram(x=ans[ans["Localization"] == organelle]["pI"], opacity=0.5, name=organelle) for organelle in organelles]
fig = {
    "data": data,
    "layout": {
        "title": "Distribution of predicted protein pI",
        "xaxis": {"title": "pH"},
        "yaxis": {"title": "Count"},
        "barmode": "overlay"
    }
}
offline.iplot(fig, show_link=False)
