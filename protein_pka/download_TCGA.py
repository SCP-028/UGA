import os
import pandas as pd


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def retrieve_annotation(df):
    """ Use manifest file to retrieve barcode information using GDC API.
    Parameters
    ----------
        df: <pandas.DataFrame>
            The manifest file used to download TCGA files.

    Returns
    -------
        annot: a <pandas.DataFrame> of more information, and
        annotDict: a <dict> of {uuid: barcode}.
    """
    url = "https://api.gdc.cancer.gov/"
    if not os.path.isfile("annot.csv"):
        if not os.path.isfile("annotation.tsv"):
            uuid = df["id"].tolist()
            url = url + "files/"
            params = {
                "filters": {
                    "op": "in",
                    "content": {
                        "field": "files.file_id",
                        "value": uuid
                    }
                },
                "format": "TSV",
                # There must be no space after comma
                "fields":
                "file_id,file_name,cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type,cases.project.project_id,cases.diagnoses.tumor_stage,cases.case_id",
                "size": len(uuid)
            }
            r = requests.post(url, json=params)  # API requires using POST method
            with open("annotation.tsv", "w") as f:
                f.write(r.text)
        annotation = pd.read_table("annotation.tsv")
        annotation = annotation[[
            "cases.0.project.project_id", "cases.0.samples.0.sample_type",
            "cases.0.diagnoses.0.tumor_stage", "cases.0.case_id",
            "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id",
            "file_id", "file_name"
        ]]
        annotation.columns = ["project", "sample_type", "tumor_stage", "caseID", "barcode", "fileID", "filename"]
        # get specific digit in barcode
        # annot.sample_type = pd.Series([int(x[-3]) for x in annot.barcode])
        # annot.loc[annot.sample_type == 0, "sample_type"] = "tumor"
        # annot.loc[annot.sample_type == 1, "sample_type"] = "normal"
        annotation.to_csv("annot.csv", index=False)
    annot = pd.read_csv("annot.csv")
    # efficiently transform to dict
    annotDict = dict(zip(annot.uuid, annot.barcode))
    return (annot, annotDict)


os.chdir("/home/jovyan/data/TCGA/count")
df = pd.read_table("/home/jovyan/data/TCGA/manifest/gdc_manifest_count.2018-08-14.txt")
url = "https://api.gdc.cancer.gov/data/"
urls = [url + ",".join(uuid) for uuid in chunks(df["id"].unique(), 100)]
with open("/home/jovyan/data/TCGA/count.urls", "w") as f:
    f.write("\n".join(urls))

"""bash
mkdir count && cd count
mv ../count.urls ./
aria2c -i count.urls
for f in ./*.gz; do tar xzf $f; done
for f in ./*/*.gz; do gunzip $f; done
cd ..
mkdir counts
"""

annot, annotDict = retrieve_annotation(df)
annot.to_csv("/home/jovyan/data/TCGA/fpkm_annot.csv", index=False)
annot = pd.read_csv("/home/jovyan/data/TCGA/annotation/fpkm_annot.csv")

annot["filename"] = annot["filename"].apply(lambda x: x[:-3] if x.endswith("gz") else x)
annot["filepath"] = annot["uuid"].str.cat(annot["filename"], sep="/")

projects = annot["project"].unique()
projects = sorted(projects)

for proj in projects:
    filenames = annot[annot["project"] == proj]["filepath"].tolist()
    ans = pd.read_table(filenames[0], header=None, index_col=0)
    for f in filenames[1:]:
        tmp = pd.read_table(f, header=None, index_col=0)
        ans = pd.concat([ans, tmp], axis=1, copy=False)
    ans.columns = [x.split("/")[0] for x in filenames]
    ans.index.names = ["Ensembl"]
    ans.to_csv(f"/home/jovyan/data/TCGA/counts/{proj}.counts.txt")
