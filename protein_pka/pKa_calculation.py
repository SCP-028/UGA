#!python3
"""
Predict protein pKa based on MCCE method.
http://pka.engr.ccny.cuny.edu/
"""
import locale
import logging
import os
import re
import subprocess
import sys
import time
from multiprocessing import Pool
from urllib.request import urlopen

import pandas as pd

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
rootpath = '/lustre1/yz73026/protein_pka'
if sys.platform == 'win32':
    os.chdir('C:/Users/jzhou/Desktop/protein_pka/')
else:
    os.chdir(rootpath)

logger = logging.getLogger('pKa_calc')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('./pdb/pKa_calculation.log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    '%(asctime)s\t%(levelname)s\t'
    '[%(filename)s:%(lineno)s -%(funcName)12s()]\t%(message)s'
    )
handler.setFormatter(formatter)
logger.addHandler(handler)


class pdb:
    def __init__(self):
        self.ids = []
        self.dl_id = []
        self.err_id = []

    def load_id(self):
        """
        Get list of existing pKa values, and list of PDB files to download
        """
        print('Loading existing pKa values...')
        annot = pd.read_csv(
            './annotation/HUMAN_9606_idmapping.dat', sep='\t', header=None)
        df = pd.read_csv('./annotation/database_charge.csv')
        idKnown = df.PDB_ID
        annot.columns = ['uniprot', 'id', 'value']
        idAll = annot.loc[annot.id == 'PDB', 'value']
        ids = list(set(idAll) - set(idKnown))
        logger.info(f'{len(ids)} PDB files need to be downloaded.')
        self.ids.append(ids)

    def getpdb(self, id, directory='pdb/'):
        """ Download PDB files from:
            ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/

        Parameters
        ----------
            id: str
                The PDB ID to download.
            directory: str, optional
                The parent directory for saving the file.

        Returns
        -------
            Nothing, just download XXXX.pdb to local directory.
        """
        pdbDir = id[1:3].lower()  # the subdirectory of the pdb files
        id = id[:4].lower()  # the pdb file names
        id_ent = f'{id}.ent'
        pdb_name = f'{id_ent}.gz'  # pdb file name
        # Make sure the download dirctory exists
        try:
            os.makedirs(os.path.join(directory, id.upper()))
        except OSError:
            pass
        saved_pdb = os.path.abspath(
            os.path.join(directory, id.upper(), f'{id.upper()}.pdb'))
        remoteaddr = f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{pdbDir}/pdb{pdb_name}'
        # print(f'Inquiring the remote file {id.upper()}.pdb ...')
        try:
            with urlopen(remoteaddr) as remotefile:
                # print(f'Saving as {saved_pdb} ...')
                with open(pdb_name, 'wb') as f:
                    f.write(remotefile.read())
            self.dl_id.append(id.upper())
            subprocess.run(['gunzip', pdb_name])
            subprocess.run(['mv', id_ent, saved_pdb])
            logger.info(f'{id.upper()} download completed.')
        except OSError:
            logger.warning(f'{id.upper()} not found.')
            self.err_id.append(id.upper())

    def preprocess(self, id, directory='pdb/', backup=True):
        """
        This program will:
        1) strip lines other than ATOM and HETATM records
        2) keep the first model of an NMR structure
        3) delete H and D atoms
        4) MSE to MET residue
        5) keep only one atom alternate position
        6) keep defined chains, if chain ID(s) are given in command
        7) remove some cofactors and salt ions

        Parameters
        ----------
            id: str
                The PDB ID to find the file.
            directory: str, optional
                The parent directory for saving the file.
            backup: bool, optional
                Whether to backup the original file or not. Default is True,
                and save to "original.bak".

        Returns
        -------
            Nothing, modify the file in place.
        """
        removable_res = [
            " ZN", "PCA", "XYP", " NA", " CL", " CA", " MG", " MN", "HOH"
        ]
        model_start = False
        newlines = []
        filepath = os.path.abspath(
            os.path.join(directory, id.upper(), f'{id.upper()}.pdb'))
        if backup:
            subprocess.run(['cp', filepath, f'{filepath}.bak'])
        with open(filepath) as f:
            for line in f:
                if line[:5] == "MODEL":
                    model_start = True
                if model_start and line[:6] == "ENDMDL":
                    break
                if line[:6] != "ATOM  " and line[:6] != "HETATM":
                    continue  # discard non ATOM records
                if line[13] == "H" or line[12] == "H":
                    continue
                if line[16] == "A":
                    line = f'{line[:16]} {line[17:]}'
                elif line[16] != " ":
                    continue  # delete this line, alternative posion is not A or empty
                if line[:6] == "HETATM" and line[17:20] == "MSE":
                    if line[12:15] == "SE ":
                        line = f"ATOM  {line[6:12]} SD{line[15:17]}MET{line[20:]}"
                    else:
                        line = f"ATOM  {line[6:17]}MET{line[20:]}"
                res = line[17:20]
                if res in removable_res:
                    continue
                newlines.append(line.rstrip())
        with open(filepath, 'w') as f:
            f.write("\n".join(newlines))
        logger.info(f'{id.upper()} preprocessing complete.')

    def set_params(self, id, directory='pdb/', quickrun=True):
        """
        Set the parameters for MCCE.

        Parameters
        ----------
            id: str
                The PDB ID of the file.
            directory: str, optional
                Directory to work in.
            quickrun: bool, optional
                Use "run.prm.quick" or "run.prm.default".

        Returns
        -------
            run.prm: a file describing the parameters that points to the PDB file.
        """
        filepath = os.path.abspath(os.path.join(directory, id.upper()))
        newlines = []
        if quickrun:
            subprocess.run([
                'cp', '/home/yz73026/src/mcce3.0/run.prm.quick',
                os.path.join(filepath, 'run.prm')
            ])
        else:
            subprocess.run([
                'cp', '/home/yz73026/src/mcce3.0/run.prm.default',
                os.path.join(filepath, 'run.prm')
            ])
        with open(os.path.join(filepath, 'run.prm')) as f:
            for line in f:
                line = line.rstrip()
                if line.endswith("(INPDB)"):
                    line = re.sub(r'^[^\s]+', fr'{id}.pdb', line)
                if line.endswith(("(DO_PREMCCE)", "(DO_ROTAMERS)",
                                  "(DO_ENERGY)", "(DO_MONTE)")):
                    line = re.sub(r'^f', r't', line)
                if line.endswith("(EPSILON_PROT)"):
                    line = re.sub(r'^[\d\.]+', r'8.0', line)
                if line.startswith("/home/mcce/mcce3.0"):
                    line = re.sub(r"^/.*3\.0", r"/home/yz73026/src/mcce3.0",
                                  line)
                newlines.append(line)
        with open(os.path.join(filepath, 'run.prm'), 'w') as f:
            f.write("\n".join(newlines))
        logger.info(f'Parameters set for {id.upper()}.')

    def calc_pka(self, id, directory='./pdb'):
        """ Calculate protein pKa values using MCCE.
            http://www.sci.ccny.cuny.edu/~jmao/mcce/manual.html

        Parameters
        ----------
            id: str
                The PDB ID of the protein calculated.
            directory: str, optional
                The parent directory for saving the file.

        Returns
        -------
            A set of files in a subdirectory named after the ID.
            See user manual for detail.
        """
        os.chdir(os.path.abspath(os.path.join(rootpath, directory, id.upper())))
        logger.info(f'{id.upper()} calculation started.')
        start = time.time()
        with open(f'{id.upper()}.run.log', 'w') as f:
            subprocess.run('/home/yz73026/src/mcce3.0/mcce', stdout=f)
        logger.info(f'{id.upper()} calculation finished, used {time.time() - start}s.')


if __name__ == '__main__':
    x = pdb()
    x.load_id()
    for item in x.ids:
        x.getpdb(item)
    subprocess.run(['find', '.', '-type', 'd', '-empty', '-delete'])
    with open('./pdb/error_pdb.list', 'w') as f:
        f.write('\n'.join(x.err_id))
    with open('./pdb/downloaded_pdb.list', 'w') as f:
        f.write('\n'.join(x.dl_id))
    for item in x.dl_id:
        x.preprocess(item)
        x.set_params(item)
    with Pool(os.cpu_count() - 1) as p:
        p.map(x.calc_pka, x.dl_id)
