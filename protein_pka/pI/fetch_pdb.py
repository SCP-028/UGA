import glob
import json
import pathlib
import subprocess
from decimal import Decimal

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline

import xmltodict as xd

# offline.init_notebook_mode()

ROOTDIR = "/home/jovyan/data/protein_pka"
pathlib.Path(f'{ROOTDIR}/result/pdb_meta').mkdir(parents=True, exist_ok=True)


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def reverse_scientic_notation(x):
    """Change string x back to PDB ID.

    Example
    -------
        '5e+28'      -> '5E28'
        '10000000.0' -> '1E07'
    """
    x = f"{Decimal(x):.0E}".split("+")
    if len(x[1]) == 1:
        x[1] = f"0{x[1]}"
    return "".join(x)


# # pK values of proteins with Ensembl IDs
if not pathlib.Path(f"{ROOTDIR}/pK_fixed.csv").is_file():
    pK = pd.read_csv(f"{ROOTDIR}/pka_cleaned_merged.csv")
    # fix IDs automatically converted to scientific notations
    pK_mask = (pK["protein"].str.len() != 4)
    pK_wrong = pK[pK_mask]
    pK.loc[pK_mask, "protein"] = pK_wrong["protein"].apply(reverse_scientic_notation)
    pK.to_csv(f"{ROOTDIR}/pK_fixed.csv", index=False)

pK = pd.read_csv(f"{ROOTDIR}/pK_fixed.csv")
print(f'There are altogether {len(pK["protein"].unique())} unique PDB IDs.')

# # Get general PDB annotation
if not pathlib.Path(f"{ROOTDIR}/result/pdb_meta/pdb_general.json").is_file():
    # downloading files using PDB's RESTful API
    url = "https://www.rcsb.org/pdb/rest/describePDB?structureId="
    pdb_urls = [url + ",".join(pdb) for pdb in chunks(pK.protein.unique(), 1000)]
    with open(f"{ROOTDIR}/result/pdb_meta/url_general.list", "w") as f:
        f.write("\n".join(pdb_urls))
    r = subprocess.run(["aria2c", "--input-file",
                        f"{ROOTDIR}/result/pdb_meta/url.list"])

    # merging downloaded files into one json
    filenames = glob.iglob(f"{ROOTDIR}/result/pdb_meta/describePDB*")
    ans = {}
    for filename in filenames:
        with open(filename, "r") as f:
            pdb = xd.parse(f.read())
        for node in pdb['PDBdescription']['PDB']:
            ans[node["@structureId"]] = node

    with open(f"{ROOTDIR}/result/pdb_meta/pdb_general.json", "w") as f:
        f.write(json.dumps(ans))

PDB_general = pd.read_json(f"{ROOTDIR}/result/pdb_meta/pdb_general.json").T
PDB_general = PDB_general[["@structureId", "@deposition_date", "@expMethod", "@resolution"]].reset_index(drop=True)
print(f"The general PDB annotation for {PDB_general.shape[0]} structures were downloaded.")

# # Get PDB entity information
if not pathlib.Path(f"{ROOTDIR}/result/pdb_meta/pdb_entity.json").is_file():
    url = "https://www.rcsb.org/pdb/rest/describeMol?structureId="
    pdb_urls = [url + ",".join(pdb) for pdb in chunks(pK.protein.unique(), 1000)]
    with open(f"{ROOTDIR}/result/pdb_meta/url_entity.list", "w") as f:
        f.write("\n".join(pdb_urls))
    r = subprocess.run(["aria2c", "--input-file",
                        f"{ROOTDIR}/result/pdb_meta/url.list"])
    filenames = glob.iglob(f"{ROOTDIR}/result/pdb_meta/describeMol*")
    ans = {}
    for filename in filenames:
        with open(filename, "r") as f:
            pdb = xd.parse(f.read())
        for node in pdb["molDescription"]["structureId"]:
            ans[node["@id"]] = node

    with open(f"{ROOTDIR}/result/pdb_meta/pdb_entity.json", "w") as f:
        f.write(json.dumps(ans))

PDB_entity = pd.read_json(f"{ROOTDIR}/result/pdb_meta/pdb_entity.json").T
print(f"The PDB entity annotation for {PDB_entity.shape[0]} structures were downloaded.")

PDB_entity_single = PDB_entity[PDB_entity["polymer"].apply(lambda x: isinstance(x, dict))]
PDB_entity_multi = PDB_entity[PDB_entity["polymer"].apply(lambda x: isinstance(x, list))]
# Unpack list in cell
s = PDB_entity_multi.apply(lambda x: pd.Series(x['polymer']), axis=1).stack().reset_index(level=1, drop=True)
s.name = "polymer"
PDB_entity_multi = PDB_entity_multi.drop("polymer", axis=1).join(s).reset_index(drop=True)
PDB_entity = pd.concat([PDB_entity_single, PDB_entity_multi], axis=0).reset_index(drop=True)
PDB_entity["@entityNr"] = PDB_entity["polymer"].apply(lambda x: x["@entityNr"])
PDB_entity["@length"] = PDB_entity["polymer"].apply(lambda x: x["@length"])
PDB_entity["@chain"] = PDB_entity["polymer"].apply(lambda x: x["chain"])  # dict
PDB_entity = PDB_entity.dropna()
print(f"After unpacking, there are {PDB_entity.shape[0]} entries in PDB_entity.")

# # Filter list, only keep human proteins
# ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/source.idx
df = pd.read_table("/home/jovyan/data/annotation/taxonomy/pdb_source.csv", header=None, sep=",")
df.columns = ["PDB", "taxonomy"]
df = df[~df["taxonomy"].isna()]
df = df[df["taxonomy"].str.contains("HOMO SAPIENS")]  # list as string

PDB_general = PDB_general[PDB_general["@structureId"].isin(df["PDB"])]
print(f"{len(PDB_general['@structureId'].unique())} entries in PDB_general are human proteins.")
PDB_entity = PDB_entity[PDB_entity["@id"].isin(df["PDB"])]
print(f"{len(PDB_entity['@id'].unique())} entries in PDB_general are human proteins.")

# PDB_general[~PDB_general["@structureId"].isin(PDB_entity["@id"].unique())]

# # Convert to Ensembl Gene ID
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
annotation = pd.read_csv("/home/jovyan/data/annotation/uniprot_id_map/HUMAN_9606_idmapping.csv")
uniprot_pdb = annotation[(annotation["ID"].isin(PDB_general["@structureId"])) & (annotation["ID_type"] == "PDB")]
uniprot_ensembl = annotation[annotation["ID_type"] == "Ensembl"]
annot = pd.merge(uniprot_pdb, uniprot_ensembl, on="UniProtKB_AC")[["ID_x", "ID_y"]].reset_index(drop=True)
annot.columns = ["PDB", "Ensembl"]
print(f'{len(pd.concat(g for _, g in annot.groupby("PDB") if len(g) > 1)["PDB"].unique())} PDB IDs map to multiple Ensembl IDs.')
print(f'{len(pd.concat(g for _, g in annot.groupby("Ensembl") if len(g) > 1)["Ensembl"].unique())} Ensembl IDs map to multiple PDB IDs.')

# # Combine previous data frames
PDB = pd.merge(PDB_entity, annot, left_on=PDB_entity["@id"], right_on=annot["PDB"])
PDB = PDB[["@id", "Ensembl", "@entityNr", "@length", "@chain"]]

PDB = pd.merge(PDB, PDB_general, left_on=PDB["@id"], right_on=PDB_general["@structureId"])
PDB = PDB.drop(["key_0", "@structureId"], axis=1).reset_index(drop=True)
print(f'There are {len(PDB["Ensembl"].unique())} unique Ensembl IDs and {len(PDB["@id"].unique())} unique PDB structures.')

# PDB[PDB["Ensembl"].duplicated()].sort_values(["Ensembl","@id", "@length", "@deposition_date", "@resolution"])
