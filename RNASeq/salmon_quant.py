#!/usr/bin/env python
# coding: utf-8

import pathlib
import subprocess

import pandas as pd

WORK_DIR = pathlib.Path("/home/yi/storage/data/JLU_mouse/")
QUANT_DIR = pathlib.Path(WORK_DIR, "results", "quant")
ANNOT_DIR = pathlib.Path("/home/yi/storage/data/annotation/genome/mouse/")

REFERENCE_GENOME_PATH = pathlib.Path(ANNOT_DIR, "GRCm38.primary_assembly.genome.fa.gz")
REFERENCE_TRANSCRIPTOME_PATH = pathlib.Path(ANNOT_DIR, "gencode.vM24.transcripts.fa.gz")
GENCODE_PATH = pathlib.Path(ANNOT_DIR, "gencode.vM24.annotation.gtf.gz")
SALMON_INDEX_DIR = pathlib.Path(ANNOT_DIR, "salmon_index")

QUANT_DIR.mkdir(parents=True, exist_ok=True)
ANNOT_DIR.mkdir(parents=True, exist_ok=True)

DOCKER_EXE = "docker exec rnaseq_Yi bash -c"

DOCKER_WORK_DIR = pathlib.Path("/home/data/data/JLU_mouse/results")
DOCKER_DATA_DIR = pathlib.Path(DOCKER_WORK_DIR, "fastp")
DOCKER_QUANT_DIR = pathlib.Path(DOCKER_WORK_DIR, "quant")
DOCKER_SALMON_IDX = pathlib.Path("/home/data/data/annotation/genome/mouse/salmon_index")


# Download files if they don't exist
if not SALMON_INDEX_DIR.exists():
    if not REFERENCE_GENOME_PATH.exists():
        subprocess.check_call(
            f"""curl -kL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz \
            -o {REFERENCE_GENOME_PATH}""",
            shell=True,
        )
    if not REFERENCE_TRANSCRIPTOME_PATH.exists():
        subprocess.check_call(
            f"""curl -kL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz \
            -o {REFERENCE_TRANSCRIPTOME_PATH}""",
            shell=True,
        )
    if not GENCODE_PATH.exists():
        subprocess.check_call(
            f"""curl -kL ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz \
            -o {GENCODE_PATH}""",
            shell=True,
        )
    subprocess.check_call(
        f'''grep "^>" <(gunzip -c {REFERENCE_GENOME_PATH}) | cut -d " " -f 1 > {ANNOT_DIR}/decoys.txt && \
            sed -i.bak -e 's/>//g' {ANNOT_DIR}/decoys.txt && \
            cat {REFERENCE_TRANSCRIPTOME_PATH} {REFERENCE_GENOME_PATH} > {ANNOT_DIR}/gentrome.fa.gz && \
            {DOCKER_EXE} "salmon index -t {ANNOT_DIR}/gentrome.fa.gz \\
            -d {ANNOT_DIR}/decoys.txt \\
            -i {ANNOT_DIR}/salmon_index \\
            --gencode -p 12"''',
        shell=True,
    )

design_mat = pd.read_csv(f"{WORK_DIR}/annotation/design_matrix.csv")
sampleIDs = design_mat.Sample.tolist()

for sID in sampleIDs:
    subprocess.check_call(
        f'''{DOCKER_EXE} "salmon quant -i {DOCKER_SALMON_IDX} -l A \
        -1 {DOCKER_DATA_DIR}/{sID}_1.fq.gz \
        -2 {DOCKER_DATA_DIR}/{sID}_2.fq.gz \
        --validateMappings --gcBias -p 8 \
        -o {DOCKER_QUANT_DIR}/{sID}_salmon"''',
        shell=True,
    )
