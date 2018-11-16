#!python3
"""
Download methylation data from the GDC (TCGA) database. Only works on Linux machines
because `uvloop` is used for asynchronous downloading.
"""
import asyncio
import glob
import logging
import os
import re
import sys

import pandas as pd
import requests

import aiohttp
import uvloop

ROOTPATH = os.path.dirname(os.path.realpath(sys.argv[0]))
os.chdir(ROOTPATH)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler(f"{__file__}.log")
handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s\t%(levelname)s\t"
    "[%(filename)s:%(lineno)s -%(funcName)12s()]\t%(message)s"
)
handler.setFormatter(formatter)
logger.addHandler(handler)


class TCGADownload:
    def __init__(self):
        self.url = "https://api.gdc.cancer.gov/"
        self.chunk_size = 65536  # 64 * 1024
        self.download = []
        self.finished = []

    def load_id(self):
        """ Load UUIDs from the manifest file and transform to download URLs.

        Returns
        -------
            Pandas dataframe containing
        """
        if not os.path.isfile("manifest.tsv"):
            url = self.url + "files"
            payload = {
                "filters": {
                    "op": "in",
                    "content": {
                        "field": "files.data_category",
                        "value": "DNA Methylation"
                    }
                },
                "format": "TSV",
                "size": 20000
            }
            r = requests.post(url, json=payload)
            with open("manifest.tsv", "w") as f:
                f.write(r.text)
        df = pd.read_table("manifest.tsv")
        logging.info(f"Manifest file contains {df.filename.count()} files.")
        # exclude existing files
        self.finished = glob.glob("./**/*.methy.*", recursive=True)
        self.finished = [x.rsplit("/", 1)[-1] for x in self.finished]
        logging.info(f"{len(self.finished)} files already exist, downloading the rest...")
        df = df[~df.filename.isin(self.finished)]
        self.download = df.id.tolist()
        return df

    def download_queue(self, urls):
        """ Create a queue to download all the given urls.

        Parameters
        ----------
            urls: list
                A list of urls to download.

        Returns
        -------
            Downloaded "*.methy.tsv" files in working directory.
        """
        logger.info(f"{len(urls)} URLs to download.")
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

    async def download_session(self, sem, work_queue):
        """ Wrapper function to link the queue and the worker.

        Parameters
        ----------
            sem: asyncio.Semaphore object.
            work_queue: asyncio.Queue object.
        """
        while not work_queue.empty():
            url = await work_queue.get()
            logger.debug(f"Got url from queue: {url}")
            async with sem:
                async with aiohttp.ClientSession() as session:
                    await self.download_worker(session, url)

    async def download_worker(self, session, url):
        """ Asynchronously download the given url."""
        uuid = url.rsplit("/", 1)[-1]
        async with session.get(url) as r:
            with open(f"./download/{uuid}.methy.tsv", "wb") as f:
                while True:
                    chunk = await r.content.read(self.chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
        logger.info(f"Downloaded UUID: {uuid}")
        self.finished.append(uuid)
        asyncio.sleep(1)

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


if __name__ == '__main__':
    x = TCGADownload()
    df = x.load_id()
    try:
        os.makedirs(os.path.join(ROOTPATH, "download"))
    except OSError:
        pass
    urls = x.url + "data/"
    x.download_queue(urls)
    annot, annotDict = x.retrieve_annotation(df)
