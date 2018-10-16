import numpy as np
import pandas as pd


def FPKM_to_TPM(df):
    """FPKM expression values to TPMs.
    The sum of all TPMs in each sample are the same.
    Refer to : https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

    Parameters
    ----------
        df: pandas.DataFrame
            Each row is a gene, and each column is a sample.

    Returns
    -------
        a pandas.DataFrame with the same shape as df.
    """
    colSum = df.sum(axis=0)
    df = df.div(colSum, axis=1) * 1e+6
    return df


pK = pd.read_csv("/home/jovyan/data/protein_pka/result/PDB_annotation.csv")
annot = pd.read_csv("/home/jovyan/data/TCGA/annotation/fpkm_annot.csv")

projects = sorted(annot["project"].unique().tolist())

for project in projects:
    geneExp = pd.read_csv(f"/home/jovyan/data/TCGA/FPKM/{project}.FPKM.csv", index_col=0)
    geneExp = FPKM_to_TPM(geneExp)
    # genes in the pK data frame
    geneExp.index = geneExp.index.str.replace(r"^(.*)\..*$", r"\1")
    geneExp = geneExp[geneExp.index.isin(pK["Ensembl"])]
    # separate normal samples
    projAnnot = annot[annot["project"] == project]
    normalUUID = projAnnot[projAnnot["sample_type"].str.lower().str.contains("normal")]["uuid"].tolist()
    normalExp = geneExp[normalUUID].mean(1)
    # separate early and late stages (and uncharacterized)
