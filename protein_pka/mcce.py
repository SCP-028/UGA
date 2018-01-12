#!python3
"""
Predict protein pKa based on MCCE method.
http://pka.engr.ccny.cuny.edu/

Require MCCE 3.0 to work: https://anaconda.org/SalahSalah/mcce/files
"""
import locale
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool
from urllib.request import urlopen

import pandas as pd

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')  # Sapelo Locale is broken, quick fix
rootpath = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(rootpath)

logger = logging.getLogger('pKa_calc')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('./pKa_calculation.log')
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
        Get list of existing pKa values, and list of PDB files to download.
        """
        logger.debug('Loading existing pKa values...')
        if not os.path.exists("./annotation/uniprot_id_mapping.dat"):
            with urlopen("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz") as remotefile:
                logger.debug('Saving UniProt ID mapping data...')
                with open("./annotation/uniprot_id_mapping.dat.gz", 'wb') as f:
                    f.write(remotefile.read())
            subprocess.run(["gunzip", "./annotation/uniprot_id_mapping.dat.gz"])
        annot = pd.read_csv('./annotation/uniprot_id_mapping.dat',
                            sep='\t', header=None,
                            names=['uniprot', 'id', 'value'])
        # df = pd.read_csv('./annotation/database_charge.csv')
        # idKnown = df.PDB_ID
        idAll = annot.loc[annot.id == 'PDB', 'value']
        try:
            x = pd.read_csv(handler.baseFilename, sep="\t", header=None,
                            names=["datetime", "level", "func", "description"])
            idFinished = x.loc[x.description.str.contains(
                               r"finished", na=False), 'description']
            idFinished = idFinished.str[0:4]
            ids = list(set(idAll) - set(idFinished))
        except Exception as e:
            logger.error(e)
            # ids = list(set(idAll) - set(idKnown))
            ids = list(set(idAll))
        logger.info(f'{len(ids)} PDB files need to be downloaded.')
        self.ids = ids

    def get_link(self, ids, directory='pdb/'):
        """ Get PDB file links from:
            ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ ,
            and create a folder to store the file.

        Parameters
        ----------
            ids: list
                The PDB IDs to download.
            directory: str, optional
                The parent directory for saving the file.

        Returns
        -------
            Links to download.
        """
        if isinstance(ids, list):
            ids = [id[:4].lower() for id in ids]  # the pdb file names
            pdb_names = [f'{id}.ent.gz' for id in ids]  # pdb file name
            pdbDirs = [id[1:3].lower() for id in ids]  # the subdirectory of the pdb files
            remoteaddr = [f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{pdbDir}/pdb{pdb_name}' for pdbDir, pdb_name in zip(pdbDirs, pdb_names)]
            # Make sure the download directory exists
            for id in ids:
                try:
                    os.makedirs(os.path.join(directory, id.upper()))
                except OSError:
                    pass
        else:
            raise TypeError(f'{id} is not a string or list.')
        return remoteaddr

    def download_file(self, url):
        """
        Parameters
        ----------
            url: str
                The url for downloading.

        Returns
        -------
            Nothing, just download and unzip the file.
        """
        id = url[-11:-7]
        ID = id.upper()
        saved_pdb = os.path.join("./pdb", ID, f'{ID}.pdb')
        try:
            with urlopen(url) as remotefile:
                logger.debug(f'Saving as {saved_pdb} ...')
                with open(f"{id}.ent.gz", 'wb') as f:
                    f.write(remotefile.read())
            # self.dl_id.append(ID)  # Global variables aren't shared in multiprocessing
            subprocess.run(['gunzip', f"{id}.ent.gz"])
            shutil.move(f"{id}.ent", saved_pdb)
            logger.info(f'{ID} download completed.')
        except OSError:
            logger.warning(f'{ID} not found.')
            # self.err_id.append(ID)

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
        filepath = os.path.realpath(
            os.path.join(directory, id.upper(), f'{id.upper()}.pdb'))
        if backup:
            shutil.copy2(filepath, f'{filepath}.bak')
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
        pkgpath = os.path.join(rootpath, "mcce3.0")
        filepath = os.path.realpath(os.path.join(directory, id.upper()))
        newlines = []
        if quickrun:
            shutil.copy2(
                os.path.join(pkgpath, "run.prm.quick"),
                os.path.join(filepath, 'run.prm')
            )
        else:
            shutil.copy2([
                os.path.join(pkgpath, "run.prm.default"),
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
                    line = re.sub(r"^/.*3\.0", pkgpath,
                                  line)
                newlines.append(line)
        with open(os.path.join(filepath, 'run.prm'), 'w') as f:
            f.write("\n".join(newlines))
        logger.info(f'Parameters set for {id.upper()}.')

    def calc_pka(self, id, directory='./pdb', clean=True):
        """ Calculate protein pKa values using MCCE.
            http://www.sci.ccny.cuny.edu/~jmao/mcce/manual.html

        Parameters
        ----------
            id: str
                The PDB ID of the protein calculated.
            directory: str, optional
                The parent directory for saving the file.
            clean: bool, optional
                Only keep the PDB file, run log and pKa output.

        Returns
        -------
            A set of files in a subdirectory named after the ID.
            See user manual for detail.
        """
        id = id.upper()
        os.chdir(os.path.realpath(os.path.join(rootpath, directory, id)))
        logger.info(f'{id} calculation started.')
        start = time.time()
        with open(f'{id}.run.log', 'w') as f:
            subprocess.run('/home/yizhou/pkg/bin/mcce3.0/mcce', stdout=f)
        with open(f'{id}.run.log', 'rb') as f:
            last = f.readlines()[-1].decode().lstrip()
        if last.startswith(("Fatal", "FATAL", "WARNING", "STOP")):
            logger.warning(f'{id} calculation aborted after {time.time() - start}s, due to {last}')
        else:
            logger.info(f'{id} calculation finished, used {time.time() - start}s.')
            shutil.move("pK.out", os.path.join(rootpath, "result", f"{id}.pka"))
        if clean:
            del_list = [i for i in os.listdir() if i not in ("pK.out", f"{id}.run.log", f"{id}.pdb.bak")]
            [os.remove(item) for item in del_list]


if __name__ == '__main__':
    for folder in ["./pdb", "./annotation"]:
        try:
            os.makedirs(folder)
        except OSError:
            pass
    x = pdb()
    x.load_id()
    urls = x.get_link(x.ids)
    with Pool(10) as p:  # TODO: Change this part to multithread
        p.map(x.download_file, urls)
    subprocess.run(['find', '.', '-type', 'd', '-empty', '-delete'])
    # Workaround for multiprocessing the download and finding downloaded files
    x.dl_id = [item for item in os.listdir(os.path.join(rootpath, "pdb")) if not item.endswith("list")]
    x.err_id = list(set(x.ids) - set(x.dl_id))
    with open('./pdb/error_pdb.list', 'w') as f:
        f.write('\n'.join(x.err_id))
    with open('./pdb/downloaded_pdb.list', 'w') as f:
        f.write('\n'.join(x.dl_id))
    if not os.path.exists(os.path.join(rootpath, "mcce3.0")):
        if not os.path.exists(os.path.join(rootpath, "mcce3.0.tar.bz2")):
            with urlopen("https://anaconda.org/SalahSalah/mcce/3.0/download/linux-64/mcce-3.0-0.tar.bz2") as remotefile:
                with open("./mcce-3.0-0.tar.bz2", 'wb') as f:
                    f.write(remotefile.read())
        subprocess.run(["tar", "-xjf", "mcce-3.0-0.tar.bz2"])
        shutil.move("./info/recipe/mcce3.0", "./mcce3.0")
        shutil.rmtree(os.path.join(rootpath, "info"), ignore_errors=True)
        shutil.rmtree(os.path.join(rootpath, "bin"), ignore_errors=True)
    for item in x.dl_id:
        try:
            x.preprocess(item)
            x.set_params(item)
        except Exception as e:
            logger.error(e)
    try:
        os.makedirs("./result")
    except OSError:
        pass
    with Pool(os.cpu_count()) as p:
        p.map(x.calc_pka, x.dl_id)
