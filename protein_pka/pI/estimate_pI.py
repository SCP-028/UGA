"""Calculate protein pI based on its pK values."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline
from scipy.stats import gaussian_kde

offline.init_notebook_mode()

pK = pd.read_csv("/home/jovyan/data/protein_pka/pK_fixed.csv")
PDB = pd.read_csv("/home/jovyan/data/protein_pka/result/PDB_annotation.csv")

# amino acids to keep
acidic = ["ASP", "GLU"]
alkaline = ["ARG", "LYS", "HIS"]
pK = pK[pK["residue"].isin(acidic + alkaline)]
pK = pK[pK["protein"].isin(PDB["@id"])]
pK = pK.sort_values(["protein", "residue"])


def eq_acid(pK, pH):
    return (-1 / (1 + 10 ** (pK - pH)))


def eq_alkali(pK, pH):
    return (1 / (1 + 10 ** (pH - pK)))


ans = pd.DataFrame(index=pK.protein.unique())
for pH in np.linspace(1, 14, 200):
    pK["charge"] = np.where(pK["residue"].isin(acidic), eq_acid(pK["pK"], pH), eq_alkali(pK["pK"], pH))
    # calculate net charge for each protein
    ans[f"pH_{pH}"] = pK.groupby("protein")["charge"].sum()
pK = pK.drop("charge", axis=1).reset_index(drop=True)
ans.to_pickle("/home/jovyan/data/protein_pka/result/net_charge_unique.pkl")

ans = pd.read_pickle("/home/jovyan/data/protein_pka/result/net_charge_unique.pkl")

zero_charge = ans.abs().idxmin(axis=1)
zero_charge = pd.to_numeric(zero_charge.str[3:])
zero_charge = pd.DataFrame({"protein": zero_charge.index, "zero_charge": zero_charge})

pK = pd.merge(pK, zero_charge, how="left", on="protein")
pK["closest"] = pK["pK"] - pK["zero_charge"]
positive_min = pK[pK["closest"] >= 0].groupby("protein")["closest"].idxmin(axis=0)
positive_min = pK.iloc[positive_min]
positive_min = positive_min[["protein", "pK"]].reset_index(drop=True)

negative_max = pK[pK["closest"] < 0].groupby("protein")["closest"].idxmax(axis=0)
negative_max = pK.iloc[negative_max]
negative_max = negative_max[["protein", "pK"]].reset_index(drop=True)

ans = pd.merge(positive_min, negative_max, on="protein")
ans["pI"] = ans.mean(axis=1)

ans.to_pickle("/home/jovyan/data/protein_pka/result/pI_unique.pkl")
ans.to_csv("/home/jovyan/data/protein_pka/result/pI_unique.csv", header=True, index=False)

ans = pd.read_pickle("/home/jovyan/data/protein_pka/result/pI_unique.pkl")
ansp = pd.merge(PDB, ans, left_on="@id", right_on="protein")
ansp = ansp[["@id", "@weight", "pI"]]
ansp = ansp.drop_duplicates("@id")

x = pd.to_numeric(ansp["pI"])
y = pd.to_numeric(ansp["@weight"])

# Calculate the point density
xy = np.vstack([x, y])
z = gaussian_kde(xy)(xy)


fig, ax = plt.subplots(figsize=(20, 10), dpi=80)
ax.scatter(x, y, c=z, s=10, edgecolor='')
plt.show()

data = [go.Histogram(x=ansp["pI"])]
fig = {
    "data": data,
    "layout": {
        "title": "Distribution of predicted protein pI",
        "xaxis": {"title": "pI"},
        "yaxis": {"title": "Count"}
    }
}
offline.iplot(fig, show_link=False)
