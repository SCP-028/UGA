#!python3
"""
Predict protein pKa based on MCCE method.
http://pka.engr.ccny.cuny.edu/

Require MCCE 3.0 to work: https://anaconda.org/SalahSalah/mcce/files
"""
import asyncio
import glob
import gzip
import locale
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool
from urllib.request import urlopen

import aioftp
import pandas as pd
import uvloop

# Sapelo Locale is broken, quick fix
locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
# Set working directory
ROOTPATH = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(ROOTPATH)
# Log settings
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler(f"./pKa_calculation_{__file__}.log")
handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s\t%(levelname)s\t"
    "[%(filename)s:%(lineno)s -%(funcName)12s()]\t%(message)s"
)
handler.setFormatter(formatter)
logger.addHandler(handler)


class pdb:
    def __init__(self):
        self.all_ids = []
        self.download_ids = []  # Download -> Unzip -> Preprocess -> Calculate
        self.unzip_ids = []  # Unzip -> Preprocess -> Calculate
        self.preprocess_ids = []  # Preprocess -> Calculate
        self.ready_ids = []  # Calculate
        self.finished_ids = []  # Successfully calculated IDs
        self.error_ids = []  # Error in download, unzip, or calculation
        # IDs this script will work on (messy queue implementation)
        self.working_ids = []

    def load_id(self):
        """
        First try to get existing pKa values,
        then get the list of PDB files to download.
        """
        for folder in ["./pdb", "./annotation", "./results"]:
            try:
                os.makedirs(folder)
            except OSError:
                pass
        self.finished_ids = [id[-8:-4] for id in glob.glob("./results/*.pka")]
        logger.debug(f"{len(self.finished_ids)} finished files.")
        # Create file even at first run so that the results folder doesn't get deleted
        with open("./results/finished_ids.list", "a") as f:
            f.write("\n".join(self.finished_ids))

        self.ready_ids = list(set(
            [id[-12:-8].upper() for id in glob.glob("./pdb/*/*.pdb.bak")]) - set(self.finished_ids))
        logger.debug(f"{len(self.ready_ids)} files ready to be calculated.")

        self.preprocess_ids = list(set([id[-8:-4].upper() for id in glob.glob(
            "./pdb/*/*.pdb") if "out" not in id]) - set(self.finished_ids) - set(self.ready_ids))
        logger.debug(
            f"{len(self.preprocess_ids)} files ready to be preprocessed.")

        self.unzip_ids = [id[-11:-7].upper() for id in glob.glob("./*.ent.gz")]
        logger.debug(f"{len(self.unzip_ids)} files ready to be unzipped.")

        if not os.path.exists("./annotation/uniprot_id_mapping.dat"):
            with urlopen("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz") as remotefile:
                logger.debug(
                    "Saving UniProt ID mapping data since it doesn't exist...")
                with open("./annotation/uniprot_id_mapping.dat.gz", "wb") as f:
                    f.write(remotefile.read())
            with gzip.open(
                "./annotation/uniprot_id_mapping.dat.gz", "rb") as inFile, open(
                    "./annotation/uniprot_id_mapping.dat", "wb") as outFile:
                shutil.copyfileobj(inFile, outFile)
            os.remove("./annotation/uniprot_id_mapping.dat.gz")
        else:
            logger.debug("UniProt ID mapping data exists.")

        logger.debug("Reading all possible PDB IDs...")
        annot = pd.read_csv("./annotation/uniprot_id_mapping.dat",
                            sep="\t", header=None,
                            names=["uniprot", "id", "value"])
        self.all_ids = annot.loc[annot.id == "PDB", "value"].tolist()
        self.download_ids = list(set(self.all_ids) - set(self.unzip_ids) - set(
            self.preprocess_ids) - set(self.ready_ids) - set(self.finished_ids))
        logger.info(
            f"{len(self.download_ids)} PDB files need to be downloaded.")

    def get_link(self, ids):
        """ Get PDB file links from:
            ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ ,
            and create folders to store the files.

        Parameters
        ----------
            ids: list
                The PDB IDs to download.

        Returns
        -------
            Links to download.
        """
        if isinstance(ids, list):
            ids = [id[:4].lower() for id in ids]  # pdb file IDs
            pdb_names = [f"{id}.ent.gz" for id in ids]  # pdb filenames
            # subdirectory of the pdb files
            pdbDirs = [id[1:3].lower() for id in ids]
            remoteaddr = [
                f"ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{pdbDir}/pdb{pdb_name}" for pdbDir, pdb_name in zip(pdbDirs, pdb_names)]
        else:
            raise TypeError(f"{id} is not a string or list.")
        return remoteaddr

    def make_dirs(self, ids):
        """Make sure the download directory exists."""
        for id in ids:
            try:
                os.makedirs(os.path.join(ROOTPATH, "pdb", id.upper()))
            except OSError:
                pass

    async def download_worker(self, session, url):
        """Download the given url to working directory."""
        url = url[len("ftp://ftp.wwpdb.org"):]
        logger.debug(f"Downloading {url}")
        try:
            await session.download(url)
            self.unzip_ids.append(url[-11:-7].upper())
        except Exception as e:
            self.error_ids.append(url[-11:-7].upper())
            logger.warning(f"Error when downloading {url}: {e}")

    async def download_session(self, sem, work_queue):
        """ Get urls from the queue and pass to worker.

        Parameters
        ----------
            sem: asyncio.Semaphore object
            work_queue: asyncio.Queue object
        """
        while not work_queue.empty():
            url = await work_queue.get()
            logger.debug(f"Got url from queue: {url}")
            async with sem:
                async with aioftp.ClientSession("ftp.wwpdb.org") as session:
                    await self.download_worker(session, url)

    def download_queue(self, urls):
        """ Create a queue to download all the given urls.

        Parameters
        ----------
            urls: list
                A list of urls to download.

        Returns
        -------
            Downloaded "*.ent.gz" files in working directory.
        """
        logger.debug(f"{len(urls)} urls to download.")
        loop = uvloop.new_event_loop()
        asyncio.set_event_loop(loop)
        q = asyncio.Queue()
        sem = asyncio.Semaphore(10)
        [q.put_nowait(url) for url in urls]
        tasks = [asyncio.ensure_future(self.download_session(sem, q))
                 for _ in range(len(urls))]
        loop.run_until_complete(asyncio.gather(*tasks))
        # Zero-sleep to allow underlying connections to close
        loop.run_until_complete(asyncio.sleep(0))
        loop.close()

    def check_mcce(self):
        """Check if MCCE 3.0 exists."""
        if not os.path.exists(os.path.join(ROOTPATH, "mcce3.0")):
            if not os.path.exists(os.path.join(ROOTPATH, "mcce3.0.tar.bz2")):
                logger.debug("MCCE isn't downloaded yet. Retrieving...")
                with urlopen("https://anaconda.org/SalahSalah/mcce/3.0/download/linux-64/mcce-3.0-0.tar.bz2") as remotefile:
                    with open("./mcce-3.0-0.tar.bz2", 'wb') as f:
                        f.write(remotefile.read())
            subprocess.run(["tar", "-xjf", "mcce-3.0-0.tar.bz2"])
            shutil.move("./info/recipe/mcce3.0", "./mcce3.0")
            shutil.rmtree(os.path.join(ROOTPATH, "info"), ignore_errors=True)
            shutil.rmtree(os.path.join(ROOTPATH, "bin"), ignore_errors=True)
        else:
            logger.info("MCCE 3.0 exists, proceeding to calculation...")

    def unzip(self, id):
        """Unzip downloaded *.ent.gz file."""
        try:
            saved_pdb = os.path.join(ROOTPATH, "pdb", id, f"{id}.pdb")
            with gzip.open(f"pdb{id.lower()}.ent.gz", "rb") as inFile, open(saved_pdb, "wb") as outFile:
                shutil.copyfileobj(inFile, outFile)
            os.remove(f"pdb{id.lower()}.ent.gz")
            self.preprocess_ids.append(id)
        except Exception as e:
            self.error_ids.append(id)
            logger.warning(f"Unzip of {id} unsuccessful: {e}")

    def preprocess(self, id, backup=True):
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
        ID = id.upper()
        filepath = os.path.join(ROOTPATH, "pdb", ID, f"{ID}.pdb")
        if backup:
            shutil.copy2(filepath, f"{filepath}.bak")
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
                    line = f"{line[:16]} {line[17:]}"
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
        with open(filepath, "w") as f:
            f.write("\n".join(newlines))
        logger.debug(f"{ID} preprocessing complete.")

    def set_params(self, id, quickrun=True):
        """
        Set the parameters for MCCE.

        Parameters
        ----------
            id: str
                The PDB ID of the file.
            quickrun: bool, optional
                Use "run.prm.quick" or "run.prm.default".

        Returns
        -------
            run.prm: a file describing the parameters that points to the PDB file.
        """
        pkgpath = os.path.join(ROOTPATH, "mcce3.0")
        ID = id.upper()
        filepath = os.path.join(ROOTPATH, "pdb", ID)
        newlines = []
        if quickrun:
            shutil.copy2(
                os.path.join(pkgpath, "run.prm.quick"),
                os.path.join(filepath, "run.prm")
            )
        else:
            shutil.copy2([
                os.path.join(pkgpath, "run.prm.default"),
                os.path.join(filepath, "run.prm")
            ])
        with open(os.path.join(filepath, "run.prm")) as f:
            for line in f:
                line = line.rstrip()
                if line.endswith("(INPDB)"):
                    line = re.sub(r"^[^\s]+", fr"{id}.pdb", line)
                if line.endswith(("(DO_PREMCCE)", "(DO_ROTAMERS)",
                                  "(DO_ENERGY)", "(DO_MONTE)")):
                    line = re.sub(r"^f", r"t", line)
                if line.endswith("(EPSILON_PROT)"):
                    line = re.sub(r"^[\d\.]+", r"8.0", line)
                if line.startswith("/home/mcce/mcce3.0"):
                    line = re.sub(r"^/.*3\.0", pkgpath,
                                  line)
                newlines.append(line)
        with open(os.path.join(filepath, "run.prm"), "w") as f:
            f.write("\n".join(newlines))
        self.ready_ids.append(ID)
        logger.debug(f"Parameters set for {ID}.")

    def split_ready_ids(self, num):
        """ A naive queue implementation for multiple scripts.

        Parameters
        ----------
            num: int
                Which part of the IDs to work on.

        Returns
        -------
            A list of the actual IDs to work on, and save the lists of IDs for
            other scripts to work with if this is the first instance.
        """
        if os.path.isfile(os.path.join(ROOTPATH, "results", "working_ids.list")):
            with open(os.path.join(ROOTPATH, "results", f"working_ids.list{num}"), "r") as f:
                self.working_ids = [line.strip() for line in f]
        else:
            n = math.ceil(len(self.ready_ids) / 10)
            self.working_ids = [self.ready_ids[i:i + n]
                                for i in range(0, len(self.ready_ids), n)]
            metafile = []
            for i, ids in enumerate(self.working_ids):
                metafile.append(os.path.join(
                    ROOTPATH, "results", f"working_ids.list{i}"))
                with open(os.path.join(ROOTPATH, "results", f"working_ids.list{i}"), "w") as f:
                    f.write("\n".join(ids))
                logger.debug(
                    f"Saved {len(ids)} IDs to file working_ids.list{i} .")
            with open(os.path.join(ROOTPATH, "results", "working_ids.list"), "w") as f:
                f.write("\n".join(metafile))
            self.working_ids = self.working_ids[num]

    def calc_pka(self, id, clean=True):
        """ Calculate protein pKa values using MCCE.
            http://www.sci.ccny.cuny.edu/~jmao/mcce/manual.html

        Parameters
        ----------
            id: str
                The PDB ID of the protein calculated.
            clean: bool, optional
                Only keep the PDB file, run log and pKa output.

        Returns
        -------
            A set of files in a subdirectory named after the ID.
            See user manual for detail.
        """
        id = id.upper()
        os.chdir(os.path.realpath(os.path.join(ROOTPATH, "pdb", id)))
        logger.info(f"{id} calculation started.")
        start = time.time()
        with open(f"{id}.run.log", "w") as f:
            subprocess.run(f"{ROOTPATH}/mcce3.0/mcce", stdout=f)
        with open(f"{id}.run.log", "rb") as f:
            last = f.readlines()[-1].decode().lstrip()
        if last.startswith(("Fatal", "FATAL", "WARNING", "STOP")):
            self.error_ids.append(id)
            logger.warning(
                f"{id} calculation aborted after {time.time() - start}s, due to {last}")
        else:
            self.finished_ids.append(id)
            logger.info(
                f"{id} calculation finished, used {time.time() - start}s.")
            shutil.move("pK.out", os.path.join(
                ROOTPATH, "results", f"{id}.pka"))
        if clean:
            del_list = [i for i in os.listdir() if i not in (
                "pK.out", f"{id}.run.log", f"{id}.pdb.bak")]
            [os.remove(item) for item in del_list]


if __name__ == "__main__":
    x = pdb()
    x.load_id()
    urls = x.get_link(x.download_ids)
    x.make_dirs(x.all_ids)
    x.download_queue(urls)

    x.check_mcce()
    for id in x.unzip_ids:
        x.unzip(id)
    for id in x.preprocess_ids:
        try:
            x.preprocess(id)
            x.set_params(id)
        except Exception as e:
            x.error_ids.append(id)
            logger.warning(f"Preprocess of {id}: {e}")
    # subprocess.run(["find", ".", "-type", "d", "-empty", "-delete"])

    x.split_ready_ids(0)  # 0 - 9, run 0 first to generate other lists

    with Pool(os.cpu_count()) as p:
        p.map(x.calc_pka, x.working_ids)

    with open("./results/finished_ids.list", "a") as f:
        f.write("\n".join(x.working_ids))

    with open("./results/error_ids.list", "a") as f:
        f.write("\n".join(x.error_ids))
