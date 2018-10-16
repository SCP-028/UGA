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


def mean_expression(df, annotation, stage_regex=r".*", normal=False):
    """Return the mean expression values of the subset of samples.

    Parameters
    ----------
        df: pandas.DataFrame
            Each row is a gene and each column is a sample.
        annotation: pandas.DataFrame
            A data frame with matching uuids and tumor stages.
        stage_regex: raw str, optional
            The regex used for filtering the "tumor_stage" column. Choose all by default.
        normal: bool, optional
            Return normal samples. Default is false.
    Returns
    -------
        a pandas.Series with names as Ensembl IDs and values as TPMs, and
        a list of sample names.
    """
    normalUUID = annotation[annotation["sample_type"].str.lower().str.contains("normal")]["uuid"].tolist()
    if normal:
        return df[normalUUID].mean(1)
    ids = annotation[
        (~annotation["uuid"].isin(normalUUID)) &
        (annotation["tumor_stage"].str.contains(stage_regex, regex=True))]["uuid"].tolist()
    return (df[ids].mean(1), ids)


def eq_acid(pK, pH):
    return (-1 / (1 + 10 ** (pK - pH)))


def eq_alkali(pK, pH):
    return (1 / (1 + 10 ** (pH - pK)))


PDB = pd.read_csv("/home/jovyan/data/protein_pka/result/PDB_annotation.csv")
pK = pd.read_csv("/home/jovyan/data/protein_pka/pK_fixed.csv")
annot = pd.read_csv("/home/jovyan/data/TCGA/annotation/fpkm_annot.csv")

acidic = ["ASP", "GLU"]
alkaline = ["ARG", "LYS", "HIS"]
pK = pK[(pK["residue"].isin(acidic + alkaline)) & pK["protein"].isin(PDB["@id"])].reset_index(drop=True)

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
    stage1Exp, stage1ID = mean_expression(df=geneExp, annotation=projAnnot, stage_regex=r"^i$|\si[abc]?$|1")
    stage2Exp, stage2ID = mean_expression(df=geneExp, annotation=projAnnot, stage_regex=r"^ii$|\si{2}[abc]?$|2")
    stage3Exp, stage3ID = mean_expression(df=geneExp, annotation=projAnnot, stage_regex=r"^iii$|\si{3}[abc]?$|3")
    stage4Exp, stage4ID = mean_expression(df=geneExp, annotation=projAnnot, stage_regex=r"^iv$|\siv[abc]?$|4")
    otherID = projAnnot[~projAnnot["uuid"].isin(normalID + stage1ID + stage2ID + stage3ID + stage4ID)]["uuid"].tolist()
    otherExp = geneExp[otherID].mean(1)
    print(f"""In project {project}, there are {projAnnot.shape[0]} samples in total:
        {len(normalID)} normal samples, {len(stage1ID)} stage i samples,
        {len(stage2ID)} stage ii samples, {len(stage3ID)} stage iii samples,
        {len(stage4ID)} stage iv samples, and {len(otherID)} samples undetermined.
    """)
