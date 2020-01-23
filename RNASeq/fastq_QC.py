#!/usr/bin/env python
# coding: utf-8

import glob
import re
import pathlib
import subprocess

import numpy as np
import pandas as pd

WORK_DIR = pathlib.Path("/home/yi/storage/data/JLU_mouse/")
DATA_DIR = pathlib.Path(WORK_DIR, "raw")
RES_DIR = pathlib.Path(WORK_DIR, "results")

DOCKER_EXE = "docker exec rnaseq_Yi bash -c"

DOCKER_WORK_DIR = pathlib.Path("/home/data/data/JLU_mouse/")
DOCKER_DATA_DIR = pathlib.Path(DOCKER_WORK_DIR, "raw")
DOCKER_RES_DIR = pathlib.Path(DOCKER_WORK_DIR, "results")

for d in [
    f"{RES_DIR}/fastp/fastp_reports",
    f"{RES_DIR}/fastp/fastqc_reports",
    f"{RES_DIR}/fastp/multiqc",
]:
    pathlib.Path(d).mkdir(parents=True, exist_ok=True)


# Generate design matrix and combine files w/ multiple lanes
all_files = glob.glob(f"{DATA_DIR}/**/FCH*.fq.gz", recursive=True)
uniq_files = set([re.sub(r"(.*)_[12]\.fq\.gz$", r"\1", x) for x in all_files])

design_mat = pd.DataFrame(
    {
        "Sample": [x.split("/")[-2] for x in uniq_files],
        "Filename": [x.split("/")[-1] for x in uniq_files],
    }
)
design_mat["SampleType"] = np.select(
    [
        design_mat["Sample"].str.startswith("M"),
        design_mat["Sample"].str.startswith("C"),
    ],
    ["Model", "Control"],
    default="Genpin",
)

design_mat = design_mat.groupby(["SampleType", "Sample"], as_index=False).agg(
    {"Filename": lambda x: ",".join(x)}
)

# e.g. a_1.fq.gz, b_1.fq.gz, c_1.fq.gz  ---> sampleA_merged_1.fq.gz
#      a_2.fq.gz, b_2.fq.gz, c_2.fq.gz  ---> sampleA_merged_2.fq.gz
sm_multi_lanes = design_mat[design_mat["Filename"].str.contains(",")]
filenames = sm_multi_lanes.apply(
    lambda x: [f"{DATA_DIR}/{x['Sample']}/{ele}" for ele in x["Filename"].split(",")],
    axis=1,
).tolist()

for f in filenames:
    _sample = pathlib.Path(f[0]).parent
    _sample_group = _sample.name
    if not pathlib.Path(_sample, f"{_sample_group}_merged_1.fq.gz").exists():
        _f_of_sample = " ".join([f"{x}_1.fq.gz" for x in f])
        _cmd = f"cat {_f_of_sample} > {_sample}/{_sample_group}_merged_1.fq.gz"
        subprocess.check_call(_cmd, shell=True)

    if not pathlib.Path(_sample, f"{_sample_group}_merged_2.fq.gz").exists():
        _f_of_sample = " ".join([f"{x}_2.fq.gz" for x in f])
        _cmd = f"cat {_f_of_sample} > {_sample}/{_sample_group}_merged_2.fq.gz"
        subprocess.check_call(_cmd, shell=True)

    design_mat.loc[
        design_mat["Sample"] == _sample_group, "Filename"
    ] = f"{_sample_group}_merged"

design_mat.to_csv(f"{WORK_DIR}/annotation/design_matrix.csv", index=False)


# Run fastp on all fastq files
filenames = [
    pathlib.Path(DOCKER_DATA_DIR, _sample, _f)
    for _sample, _f in zip(design_mat.Sample, design_mat.Filename)
]
for f, sg in zip(filenames, design_mat.Sample):
    if pathlib.Path(RES_DIR, "fastp", "fastp_reports", f"{sg}_report.json").exists():
        continue
    subprocess.check_call(
        f'''{DOCKER_EXE} "fastp -V -i {f}_1.fq.gz -I {f}_2.fq.gz \
            -o {DOCKER_RES_DIR}/fastp/{sg}_1.fq.gz \
            -O {DOCKER_RES_DIR}/fastp/{sg}_2.fq.gz \
            --html {DOCKER_RES_DIR}/fastp/fastp_reports/{sg}_fastp.html \
            --json {DOCKER_RES_DIR}/fastp/fastp_reports/{sg}_fastp.json \
            -w 4"''',
        shell=True,
    )

subprocess.check_call(
    f'''{DOCKER_EXE} "multiqc {DOCKER_RES_DIR}/fastp/fastp_reports/ \
        -m fastp -o {DOCKER_RES_DIR}/fastp/multiqc/fastp/"''',
    shell=True,
)


# ## Run FastQC on all fastp-generated files
if not pathlib.Path(
    RES_DIR, "fastp", "multiqc", "fastqc", "multiqc_report.html"
).exists():
    subprocess.check_call(
        f'''{DOCKER_EXE} "fastqc {DOCKER_RES_DIR}/fastp/*.fq.gz \
            --noextract -o {DOCKER_RES_DIR}/fastp/fastqc_reports/ -t 8"''',
        shell=True,
    )
    subprocess.check_call(
        f'''{DOCKER_EXE} "multiqc {DOCKER_RES_DIR}/fastp/fastqc_reports/ \
            -m fastqc -o {DOCKER_RES_DIR}/fastp/fastqc_reports/"''',
        shell=True,
    )
