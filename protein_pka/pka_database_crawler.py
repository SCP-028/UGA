#!python3
"""
http://pka.engr.ccny.cuny.edu/index.php
Crawler for the pKa Database, retrieving all the pKa values for the calculated PDB proteins.
"""
import re
from io import StringIO
import requests
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup


def ParseTable(soup, id):
    """Format the retrieved text table into a DataFrame.

    Parameters
    ----------
        soup: BeautifulSoup
            The HTML table to be formatted.
        id: str
            The PDB ID, used for cleaning possible mistakes in the database.

    Return
    ------
        A pandas.DataFrame consisting of four columns:
        Record_ID, PDB_ID, Residue, pKa
    """
    df = soup.find('pre').string
    df = re.sub(r', ', r',', df)  # Probably multiple chains in Record_ID
    df = re.sub(r' +', r'\t', df)
    df = re.sub(r'\n\t', r'\n', df)  # multiple spaces -> tabs
    df = re.sub(r'[^\t]*=.*', 'NA', df)  # Pairwise values removed
    pdb = re.compile(f'({id})_?([^\t\.]+)')
    df = re.sub(pdb, r'\1\t\2', df)  # fix PDB_IDs that have no spaces after
    df = re.sub(r'<|>', r'', df)
    df = StringIO(df)
    # df.seek(0)  # Seek back to the start of the file
    datatype = {
        'Record_ID': 'str',
        'PDB_ID': 'str',
        'Residue': 'str',
        'pKa': np.float32
    }
    df = pd.read_table(
        df,
        sep='\t',
        skiprows=1,
        usecols=['Record_ID', 'PDB_ID', 'Residue', 'pKa'],
        dtype=datatype)
    return df


# Get all PDB IDs in the database
HEADER = {
    'User-Agent':
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/62.0.3202.75 Safari/537.36',
    'Accept':
    'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
    'DNT':
    '1',
    'Referer':
    'http://pka.engr.ccny.cuny.edu/index.php',
    'Accept-Encoding':
    'gzip, deflate',
    'Accept-Language':
    'zh-CN,zh;q=0.9,en;q=0.8,zh-TW;q=0.7'
}
URL = 'http://pka.engr.ccny.cuny.edu/browse.php'
SS = requests.Session()
r = SS.get(URL, headers=HEADER)
soup = BeautifulSoup(r.content, 'lxml')
PROTEINS = soup.find_all('a', href=re.compile(r'browse.php\?protein=.+'))
PROTEINS = [x.string for x in PROTEINS]
PDB = []
regex = re.compile(r'browse.php\?protein=.+pdb=\w{4}$')
for protein in PROTEINS:
    r = SS.get(URL, headers=HEADER, params={'protein': protein})
    soup = BeautifulSoup(r.content, 'lxml')
    PDB.append([x.string for x in soup.find_all('a', href=regex)])
PDBIDS = [item for sublist in PDB for item in sublist]

# Retrieve pKa values using all the PDB IDs
URL = 'http://pka.engr.ccny.cuny.edu/get_txt_result.php'
del HEADER['Referer']
PARAM = {
    'residue': 'ALL',
    'name_type': '1',
    'Submit': 'Go',
    'method': 'ALL',
    'pkamethod': 'ALL',
    'range_min': '0',
    'range_max': '14',
    'resolution': '3.0',
    'size_min': '0',
    'size_max': '99999',
    'dsolv_min': '0',
    'dsolv_max': '99.00',
    'submitter': ''
}
pKa = pd.DataFrame(columns=['Record_ID', 'PDB_ID', 'Residue', 'pKa'])
for id in PDBIDS:
    PARAM['protein_string'] = id
    r = SS.get(URL, headers=HEADER, params=PARAM)
    soup = BeautifulSoup(r.content, 'lxml')
    df = ParseTable(soup, id)
    pKa = pKa.append(df)

pKa.to_csv("./pKa.csv", index_label=False, index=False)
