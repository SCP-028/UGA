"""Calculate protein pI based on its pK values."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline

offline.init_notebook_mode()

df = pd.read_csv("/data/protein_pka/pka_cleaned_merged.csv")
with open("/data/protein_pka/pdb_current.json", "r") as f:
    pdb = json.load(f)

# amino acids to keep
acidic = ["ASP", "GLU"]
alkaline = ["ARG", "LYS", "HIS"]
df = df[df["residue"].isin(acidic + alkaline)]
df = df[df["protein"].isin(pdb["idList"])]
df = df.sort_values(["protein", "residue"])


def eq_acid(pK, pH):
    return (-1 / (1 + 10 ** (pK - pH)))


def eq_alkali(pK, pH):
    return (1 / (1 + 10 ** (pH - pK)))


ans = pd.DataFrame(index=df.protein.unique())
for pH in np.linspace(1, 14, 100):
    df["charge"] = np.where(df["residue"].isin(acidic), eq_acid(df["pK"], pH), eq_alkali(df["pK"], pH))
    # calculate net charge for each protein
    ans[f"pH_{pH}"] = df.groupby("protein")["charge"].sum()
df = df.drop("charge", axis=1)
df = df.reset_index(drop=True)

ans.to_pickle("/data/protein_pka/result/net_charge.pkl")

zero_charge = ans.abs().idxmin(axis=1)
zero_charge = pd.to_numeric(zero_charge.str[3:])
zero_charge = pd.DataFrame({
    "protein": zero_charge.index, "zero_charge": zero_charge})

data = [go.Histogram(x=zero_charge["zero_charge"])]
fig = {
    "data": data,
    "layout": {
        "title": "Distribution of zero-charged proteins' pH level",
        "xaxis": {"title": "pH"},
        "yaxis": {"title": "Count"}
    }
}
offline.iplot(fig, show_link=False)

df = pd.merge(df, zero_charge, how="left", on="protein")

df["closest"] = df["pK"] - df["zero_charge"]
positive_min = df[df["closest"] >= 0].groupby("protein")["closest"].idxmin(axis=0)
positive_min = df.iloc[positive_min]
positive_min = positive_min[["protein", "pK"]].reset_index(drop=True)

negative_max = df[df["closest"] < 0].groupby("protein")["closest"].idxmax(axis=0)
negative_max = df.iloc[negative_max]
negative_max = negative_max[["protein", "pK"]].reset_index(drop=True)

ans = pd.merge(positive_min, negative_max, on="protein")
ans["pI"] = ans.mean(axis=1)
ans.to_pickle("/data/protein_pka/result/pI.pkl")

data = [go.Histogram(x=ans["pI"])]
fig = {
    "data": data,
    "layout": {
        "title": "Distribution of predicted protein pI",
        "xaxis": {"title": "pH"},
        "yaxis": {"title": "Count"}
    }
}
offline.iplot(fig, show_link=False)
