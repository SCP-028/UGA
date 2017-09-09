#!python3
import os
import errno
import re
import glob
import gzip
import shutil
import requests
import sqlite3
import pandas as pd

os.chdir("C:/Users/jzhou/Desktop/raw_expression")


def downloadData(df, directory='./sep'):
    """Use manifest file to download data using GDC data API.
    df: [pandas data frame] of the manifest file downloaded from GDC website.
    """
    homeDir = os.getcwd()
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    os.chdir(directory)
    uuid = df.id.tolist()
    url = 'https://api.gdc.cancer.gov/data/'
    uuid = [url + x for x in uuid]
    fileNum = len(uuid)
    for id in uuid:
        os.system(f"curl --remote-name --remote-header-name '{id}")
    print(f"Downloaded {fileNum} files to './sep/ .")
    os.chdir(homeDir)


def uuidToBarcode(df):
    """Use manifest file to retrieve barcode information using GDC API.
    df: a [pandas dataframe] of the manifest file used to download TCGA files.
    Returning annot: a [pandas dataframe] of more information, and
    annotDict: a dict of {filename: barcode}.
    """
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
        "fields": "file_id,file_name,cases.samples.submitter_id,cases.samples.sample_type,cases.project.project_id,cases.diagnoses.tumor_stage",
        "size": len(uuid)
    }
    url = "https://gdc-api.nci.nih.gov/files"
    r = requests.post(url, json=params)  # API requires using POST method
    with open("./results/annotation.tsv", "w") as f:
        f.write(r.text)  # save raw annotation file
    annot = pd.read_csv("./results/annotation.tsv", sep="\t")
    annot = annot[['file_name', 'cases.0.project.project_id',
                   'cases.0.samples.0.submitter_id',
                   'cases.0.samples.0.sample_type']]
    annot = annot.rename(columns={
        'cases.0.project.project_id': 'project',
        'cases.0.samples.0.submitter_id': 'barcode',
        'cases.0.samples.0.sample_type': 'sample_type'
    })
    annot.file_name = annot.file_name.str.replace(
        '.gz', '')  # regex in pandas dataframe
    annot.project = annot.project.str.replace('TCGA.', '')
    annot.sample_type = annot.sample_type.str.replace(
        '^.*normal.*$', 'normal', case=False)  # set all values containing 'normal' to normal, and everything else to tumor
    annot.sample_type[~annot.sample_type.isin(['normal'])] = 'tumor'
    # efficiently transform to dict
    annotDict = dict(zip(annot.file_name, annot.barcode))

    return(annot, annotDict)


def unzipAll():
    """Unzip all txt.gz files downloaded by the GDC file transfer tool.
    WARNING: will remove all zipfiles!
    """
    for zipfile in glob.iglob('./**/*.gz', recursive=True):
        newfile = re.sub('.gz$', '', zipfile)
        with gzip.open(zipfile, 'rb') as f_in, open(newfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(zipfile)


def mergeData(annot, annotDict):
    """Merge all the downloaded data by column
    annot: [pandas dataframe] annot from function {uuidToBarcode}
    annotDict: [dict] from function {uuidToBarcode}
    return merged [pandas dataframe] df
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
                './sep/' + cases[0], sep='\t', names=['ensembl', annotDict[cases[0]]])
            cases.pop(0)  # Get first case (ensembls) and remove it from list
            geneNum = len(df)
            for case in cases:
                dfSingle = pd.read_csv(
                    './sep/' + case, sep="\t", names=['ensembl', annotDict[case]])
                # df = pd.merge(df, dfSingle, how='outer', on='ensembl')
                if geneNum == len(dfSingle):
                    df = pd.concat([df, dfSingle.iloc[:, 1]], axis=1)
                else:
                    df = pd.merge(df, dfSingle, how='outer', on='ensembl')
            # df.to_sql(name=project + '_normal', con=db, if_exists='replace')
            df.to_csv('./results/' + project +
                      '_normal.csv', sep='\t', index=False)
        print(project + ' normal finished!')
        # Tumor
        annotation = annot[(annot.project == project) &
                           (annot.sample_type == "tumor")]
        cases = annotation.file_name.tolist()
        if len(cases) != 0:
            df = pd.read_csv(
                './sep/' + cases[0], sep='\t', names=['ensembl', annotDict[cases[0]]])
            cases.pop(0)
            geneNum = len(df)
            for case in cases:
                dfSingle = pd.read_csv(
                    './sep/' + case, sep="\t", names=['ensembl', annotDict[case]])
                if geneNum == len(dfSingle):
                    df = pd.concat([df, dfSingle.iloc[:, 1]], axis=1)
                else:
                    df = pd.merge(df, dfSingle, how='outer', on='ensembl')
            # df.to_sql(name=project + '_tumor', con=db, if_exists='replace')
            df.to_csv('./results/' + project +
                      '_tumor.csv', sep='\t', index=False)
            print(project + ' tumor finished!')
    # db.close()


def run(manifest="gdc_manifest.2017-08-22T20-13-51.839644.txt"):
    df = pd.read_csv(manifest, sep="\t")
    # downloadData(df)
    annot, annotDict = uuidToBarcode(df)
    # unzipAll()
    mergeData(annot, annotDict)


if __name__ == '__main__':
    run()
