#!python3
"""
Download gene expression data from the GDC (TCGA) database.
"""
import os
import errno
import logging
import re
import glob
import gzip
import shutil
import requests
import pandas as pd

logging.basicConfig(filename='./annotation/download.log', level=logging.INFO)
try:
    os.chdir("/home/yizhou/dockers/RStudio/data/expression_count")
except BaseException:
    os.chdir("C:/users/jzhou/Desktop/expression_count")


def downloadData(df, directory='./sep'):
    """Use manifest file to download data using GDC data API.
    Arguements
        df: [pandas data frame] of the manifest file downloaded from GDC website.
        directory: a [str] showing the directory to store the downloaded data
    """
    homeDir = os.getcwd()
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    os.chdir(directory)
    fileNum = df.filename.count()
    logging.info(f"Manifest file contains {fileNum} files.")
    # exclude existing files
    # change counts to FPKM if downloading FPKM data
    existFile = glob.glob("./**/*.counts.*", recursive=True)
    existFile = [
        re.sub(r".*\/(.*\.txt)(\.gz)?$", r"\1.gz", x) for x in existFile
    ]  # include unzipped files
    fileNum = len(existFile)
    logging.info(f"{fileNum} files already exist, downloading the rest...")
    url = 'https://api.gdc.cancer.gov/data/'
    df = df[~df.filename.isin(existFile)]
    # download files
    uuid = df.id.tolist()
    uuid = [url + x for x in uuid]
    fileNum = len(uuid)
    for id in uuid:
        os.system(f"curl --remote-name --remote-header-name {id}")
    logging.info(f"Downloaded {fileNum} files to {directory}")
    os.chdir(homeDir)


def uuidToBarcode(df, directory='./annotation'):
    """Use manifest file to retrieve barcode information using GDC API.
    Arguments
        df: a [pandas dataframe] of the manifest file used to download TCGA files.
        directory: a [str] showing the directory to store annotation.tsv and annot.tsv

    Return
        annot: a [pandas dataframe] of more information, and
        annotDict: a dict of {filename: barcode}.
    """
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    annotFile = glob.glob(f"{directory}/annotation.tsv", recursive=True)
    if not annotFile:
        uuid = df.id.tolist()
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
            "fields": "file_id,file_name,cases.samples.submitter_id,cases.samples.sample_type,cases.project.project_id,cases.diagnoses.tumor_stage,cases.case_id",
            "size": len(uuid)
        }
        url = "https://api.gdc.cancer.gov/files"
        r = requests.post(url, json=params)  # API requires using POST method
        with open(f"{directory}/annotation.tsv", "w") as f:
            f.write(r.text)  # save raw annotation file
    annot = pd.read_table(f"{directory}/annotation.tsv")
    annot = annot[[
        'file_name', 'cases.0.project.project_id',
        'cases.0.samples.0.submitter_id', 'cases.0.samples.0.sample_type',
        'cases.0.diagnoses.0.tumor_stage'
    ]]
    annot = annot.rename(columns={
        'cases.0.project.project_id': 'project',
        'cases.0.samples.0.submitter_id': 'barcode',
        'cases.0.samples.0.sample_type': 'sample_type',
        'cases.0.diagnoses.0.tumor_stage': 'tumor_stage'
    })
    annot.file_name = annot.file_name.str.replace(
        '.gz', '')  # regex in pandas dataframe
    annot.project = annot.project.str.replace('TCGA.', '')
    # get specific digit in barcode
    annot.sample_type = pd.Series([int(x[-3]) for x in annot.barcode])
    annot.loc[annot.sample_type == 0, 'sample_type'] = 'tumor'
    annot.loc[annot.sample_type == 1, 'sample_type'] = 'normal'
    annot.to_csv(f"{directory}/annot.tsv", index=False)
    # efficiently transform to dict
    annotDict = dict(zip(annot.file_name, annot.barcode))

    return (annot, annotDict)


def unzipAll():
    """Unzip all txt.gz files downloaded by the GDC file transfer tool.
    WARNING: will remove all zipfiles!
    """
    for zipfile in glob.iglob('./**/*.gz', recursive=True):
        newfile = re.sub('.gz$', '', zipfile)
        with gzip.open(zipfile, 'rb') as f_in, open(newfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(zipfile)


def mergeData(annot, annotDict, directory="./results", filedir='./sep'):
    """Merge all the downloaded data by column
    Arguement
        annot: [pandas dataframe] annot from function `uuidToBarcode`
        annotDict: [dict] from function `uuidToBarcode`
    """
    # db = sqlite3.connect('./results/results.sql')
    projects = annot.project.unique().tolist()
    for project in projects:
        # Normal
        annotation = annot[(annot.project == project) &
                           (annot.sample_type == "normal")]
        cases = annotation.file_name.tolist()
        if len(cases) != 0:
            df = pd.read_csv(
                f'{filedir}/{cases[0]}',
                sep='\t',
                names=['ensembl', annotDict[cases[0]]])
            cases.pop(0)  # Get first case (ensembls) and remove it from list
            for case in cases:
                try:
                    dfSingle = pd.read_csv(
                        f'{filedir}/{case}',
                        sep="\t",
                        names=['ensembl', annotDict[case]])
                    df = pd.merge(df, dfSingle, how='outer', on='ensembl')
                except FileNotFoundError as e:
                    logging.warning(e)
            # df.to_sql(name=project + '_normal', con=db, if_exists='replace')
            df.to_csv(
                f"{directory}/{project}_normal.csv", sep='\t', index=False)
        logging.info(f'{project} normal finished!')
        # Tumor
        annotation = annot[(annot.project == project)
                           & (annot.sample_type == "tumor")]
        cases = annotation.file_name.tolist()
        if len(cases) != 0:
            df = pd.read_csv(
                f'{filedir}/{cases[0]}',
                sep='\t',
                names=['ensembl', annotDict[cases[0]]])
            cases.pop(0)
            for case in cases:
                try:
                    dfSingle = pd.read_csv(
                        f'{filedir}/{case}',
                        sep="\t",
                        names=['ensembl', annotDict[case]])
                    df = pd.merge(df, dfSingle, how='outer', on='ensembl')
                except FileNotFoundError as e:
                    logging.warning(e)
            # df.to_sql(name=project + '_tumor', con=db, if_exists='replace')
            df.to_csv(
                f"{directory}/{project}_tumor.csv", sep='\t', index=False)
            logging.info(f'{project} tumor finished!')
    # db.close()


def run(manifest="manifest.txt"):
    df = pd.read_csv(manifest, sep="\t")
    downloadData(df)
    annot, annotDict = uuidToBarcode(df)
    unzipAll()
    mergeData(annot, annotDict)


if __name__ == '__main__':
    run()
