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


def unpack_list_in_df(df, col_target):
    # Flatten columns of lists
    col_flat = [item for sublist in df[col_target] for item in sublist]
    # Row numbers to repeat
    lens = df[col_target].apply(len)
    vals = range(df.shape[0])
    ilocations = np.repeat(vals, lens)
    # Replicate rows and add flattened column of lists
    cols = [i for i, c in enumerate(df.columns) if c != col_target]
    new_df = df.iloc[ilocations, cols].copy()
    new_df[col_target] = col_flat
    return new_df


# # pK values of proteins with Uniprot IDs
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
PDB_general = PDB_general[PDB_general["@status"] == "CURRENT"]
PDB_general = PDB_general[["@structureId", "@deposition_date", "@expMethod", "@resolution"]].reset_index(drop=True)
print(f"The general PDB annotation for {PDB_general.shape[0]} structures were downloaded.")
# PDB_general.sample(1)

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
PDB_entity_multi = unpack_list_in_df(PDB_entity_multi, "polymer")

PDB_entity = pd.concat([PDB_entity_single, PDB_entity_multi], axis=0).reset_index(drop=True)
PDB_entity["@length"] = PDB_entity["polymer"].map(lambda x: x["@length"])
PDB_entity["@chain"] = PDB_entity["polymer"].map(lambda x: x["chain"])  # dict
PDB_entity["@weight"] = PDB_entity["polymer"].map(lambda x: x["@weight"])
PDB_entity["@Taxonomy"] = [x["Taxonomy"] if "Taxonomy" in x else np.nan for x in PDB_entity["polymer"]]
PDB_entity["@synonym"] = [x["synonym"] if "synonym" in x else np.nan for x in PDB_entity["polymer"]]
PDB_entity["@uniprot"] = [x["macroMolecule"] if "macroMolecule" in x else np.nan for x in PDB_entity["polymer"]]
PDB_entity = PDB_entity.drop("polymer", axis=1).dropna()

# represent `@chain` in a way that's easier to manipulate
PDB_entity["@chain"] = [[x] if isinstance(x, dict) else x for x in PDB_entity["@chain"]]
PDB_entity["@chain"] = PDB_entity["@chain"].map(lambda x: "".join(sorted([y["@id"] for y in x])))
PDB_entity["@chain"] = PDB_entity["@chain"].str.upper()

# same goes for `@synonym`
PDB_entity["@synonym"] = [[x] if isinstance(x, dict) else x for x in PDB_entity["@synonym"]]
PDB_entity["@synonym"] = PDB_entity["@synonym"].map(lambda x: "".join(sorted([y["@name"] for y in x])))
PDB_entity["@synonym"] = PDB_entity["@synonym"].str.upper()

# and uniprot
PDB_entity["@uniprot"] = [[x] if isinstance(x, dict) else x for x in PDB_entity["@uniprot"]]
PDB_entity["@uniprot"] = PDB_entity["@uniprot"].map(lambda x: [y["accession"]["@id"] for y in x])
PDB_entity = unpack_list_in_df(PDB_entity, "@uniprot")
print(f"After unpacking, there are {PDB_entity.shape[0]} entries in PDB_entity.")

# # Convert to Ensembl Gene ID
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
annotation = pd.read_csv("/home/jovyan/data/annotation/uniprot_id_map/HUMAN_9606_idmapping.csv")
uniprot_ensembl = annotation[annotation["ID_type"] == "Ensembl"]
PDB = pd.merge(PDB_entity, uniprot_ensembl, left_on=PDB_entity["@uniprot"], right_on=uniprot_ensembl["UniProtKB_AC"])

# # Combine previous data frames
PDB = PDB[PDB["@Taxonomy"].astype(str).str.contains("Homo sapiens")]
PDB = PDB[["@uniprot", "@id", "ID", "@length", "@weight", "@chain", "@synonym"]]
PDB = pd.merge(PDB, PDB_general, left_on=PDB["@id"], right_on=PDB_general["@structureId"])
PDB = PDB.drop(["key_0", "@structureId"], axis=1).reset_index(drop=True)
PDB = PDB.rename(columns={"ID": "Ensembl"})

print(f'There are {len(PDB["Ensembl"].unique())} unique Ensembl IDs and {len(PDB["@id"].unique())} unique PDB structures.')
print(f'{len(PDB_entity[~PDB_entity["@id"].isin(PDB["@id"])]["@id"].unique())} PDB entries don\'t have matching Ensembl IDs / are\'t human proteins.')

# # Remove duplicates
# Same gene, same structure, same chains (technical replicates)
PDB = PDB.sort_values("@length", ascending=False).drop_duplicates(keep="first")
# each unique protein, only keep the ones with max(@length)
PDB = PDB[PDB["@length"] == PDB.groupby(["@uniprot", "@chain", "Ensembl"])["@length"].transform(max)]
# remaining duplicates, keep the one with latest deposition date
PDB = PDB[PDB["@deposition_date"] == PDB.groupby(["@uniprot", "@chain", "Ensembl"])["@deposition_date"].transform(max)]
# ...highest resolution
PDB = PDB[PDB["@resolution"] == PDB.groupby(["@uniprot", "@chain", "Ensembl"])["@resolution"].transform(min)]
print(f'There are {len(PDB["@id"].unique())} unique PDB structures, and {len(PDB["Ensembl"].unique())} unique Ensembl genes')
