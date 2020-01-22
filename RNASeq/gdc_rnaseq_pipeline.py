""" RNA-Seq analysis pipeline following the methods used by ICGC.

This script works specifically for the mice data on the JLU server.
Some of the functions could be rewritten for general use.

Function:
    1. Walk through DATA_DIR and generate a design matrix.
    2. Concatenate samples with multiple fastq files (splitted lanes).
    3. Run fastp on all fq.gz files, and store the trimmed file under RES_DIR/data.
        Reports of fastp are stored in RES_DIR/fastp.
    4. Align each sample using a two-pass method with STAR

Output:
    Under `RES_DIR/data` and `RES_DIR/fastp`, we have the trimmed FASTQ files and their
    corresponding fastp reports. MultiQC reports is under `RES_DIR/fastp/multiqc`.

    Under `RES_DIR/bam`, each sample should have its own sub-directory containing
    the following files with the `sample_group` as a prefix:

    - Aligned.out.bam: all genomic alignments including chimeric and unaligned reads
    - Aligned.toTranscriptome.out.bam: aligned reads with transcript coordinates rather than genomic coordinates
    - Chimeric.out.junction: reads that were mapped to different chromosomes or strands (fusion alignments)
    - SJ.out.tab: high confidence collapsed splice junctions
    - Log(.final|.progress).out

    Under `RES_DIR/counts`, the counts produced by STAR is moved here.
    See https://www.biostars.org/p/218995/. In our case cols 1 and 2 should be kept.

    Under `RES_DIR/tpm`, TPM values produced by Salmon is stored.
    Ensembl transcript IDs are used because it's mapped to the reference transcriptome.

Software and data:
    - fastp v0.20.0: https://github.com/OpenGene/fastp
    - FastQC v0.11.9: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    - multiQC v1.8: https://multiqc.info
    - STAR v2.7.2b: https://github.com/alexdobin/STAR
    - Salmon v1.1.0: https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz
    - Reference genome: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz
    - Reference transcriptome: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz
    - Gene annotation: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz

References:
    https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
    https://github.com/akahles/icgc_rnaseq_align/blob/master/star_align.py
    https://salmon.readthedocs.io/en/latest/salmon.html
    https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    Choice of software: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/
"""
import glob
import gzip
import logging
import os
import re
import shutil
import subprocess
from typing import List

import pandas as pd

# Filepath variables
WORK_DIR = os.path.expanduser("~/rna_seq")
DATA_DIR = os.path.join(WORK_DIR, "mm_liver_data")
RES_DIR = os.path.join(WORK_DIR, "results")

REFERENCE_GENOME_PATH = os.path.join(WORK_DIR, "genome", "GRCm38.p6.genome.fa")
REFERENCE_TRANSCRIPTOME_PATH = os.path.join(
    WORK_DIR, "genome", "gencode.vM24.transcripts.fa"
)
GENCODE_PATH = os.path.join(WORK_DIR, "genome", "gencode.vM24.annotation.gtf")
STAR_INDEX_DIR = os.path.join(WORK_DIR, "star_index")

FASTP_PATH = os.path.expanduser("~/pkg/bin/fastp")
FASTQC_PATH = os.path.expanduser("~/pkg/bin/fastqc")
MULTIQC_PATH = os.path.expanduser("~/.local/bin/multiqc")
STAR_PATH = os.path.expanduser("~/pkg/bin/STAR")
SALMON_PATH = os.path.expanduser("~/pkg/bin/salmon/bin/salmon")

for d in [
    f"{RES_DIR}/data",
    f"{RES_DIR}/fastp/multiqc",
    f"{RES_DIR}/fastqc",
    f"{RES_DIR}/bam",
    f"{RES_DIR}/counts",
    f"{RES_DIR}/tpm",
    STAR_INDEX_DIR,
]:
    os.makedirs(d, exist_ok=True)

# Log settings
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler(os.path.join(f"{__file__.rstrip('.py')}.log"))
handler.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s [%(levelname)s]{%(lineno)d:%(funcName)12s()} %(message)s",
    "%Y-%m-%d %H:%M:%S",
)
handler.setFormatter(formatter)
logger.addHandler(handler)


def ListFastqFiles(group: str = "original") -> List[List[str]]:
    """ List fastq files to be analyzed in the current working directory.

    Folders could contain multiple files for the same sample,
    which should be joined as they are from multiple lanes.

    Returns:
        A list of lists of absolute paths to the fq.gz files. Each sublist
        contains files in the same directory.
    """
    # Find the lowest-level directories containing the files
    all_files = glob.glob(f"{DATA_DIR}/{group}/**/FCH*.fq.gz", recursive=True)
    dir_names = set((os.path.dirname(x) for x in all_files))
    res = []
    for dir_name in dir_names:
        # Assuming we don't have any missing files
        uniq_samples = list(
            set(
                (
                    re.sub(r"^(.*)_[12]\.fq\.gz$", r"\1", x)
                    for x in glob.glob(f"{dir_name}/FCH*.fq.gz")
                )
            )
        )
        logger.debug(f"Found {len(uniq_samples)} for {dir_name}")
        res.append(uniq_samples)
    return res


def Samples2DesignMatrix(samples: List[List[str]]) -> pd.DataFrame:
    """ Convert the samples list to a design matrix.

    Args:
        samples (List[List[str]]): return from function `ListFastqFiles`.

    Returns:
        pd.DataFrame: a design matrix with the following columns:
            - timepoint: 1 - 9
            - sample_type: Control / Model / Genpin
            - sample_group: e.g. original/T1_C/11_90bp
            - filename: comma-separated filenames (_[1|2] removed)
    """
    _samples = [x for subl in samples for x in subl]
    _samples = [re.sub(f"{DATA_DIR}/", "", x) for x in _samples]
    filenames = [os.path.basename(x) for x in _samples]
    sample_groups = [os.path.dirname(x) for x in _samples]
    possible_sample_types = {"C": "Control", "M": "Model", "G": "Genpin"}
    timepoints = []
    sample_type = []
    for sg in sample_groups:
        if sg.startswith("original"):
            tmp = re.search("/T(\d+)_([CMG])", sg)
            timepoints.append(tmp.group(1))
            sample_type.append(possible_sample_types[tmp.group(2)])

        else:  # A "new" sample
            tmp = re.search("/Time_(\d+)_.*/(Control|Model)", sg)
            timepoints.append(tmp.group(1))
            sample_type.append(tmp.group(2))
    df = pd.DataFrame(
        {
            "timepoint": timepoints,
            "sample_type": sample_type,  # Control/Model
            "sample_group": sample_groups,  #
            "filename": filenames,  # FCH*_[1|2].fq.gz
        }
    )
    df = df.sort_values(by=["timepoint", "sample_type", "sample_group"])
    df = (
        df.groupby(["timepoint", "sample_type", "sample_group"])["filename"]
        .agg(",".join)
        .reset_index()
    )
    return df


def ConcatSamples(samples: List[str]):
    """ Concatenate samples in the same group to two fq.gz files.

    Files are cat together because they are the same sample split to multiple
    lanes. The concatenate files are for the forward and reverse reads.

    Args:
        samples (List[str]): comma-separated absolute file paths with the
            trailing _[1|2].fq.gz stripped

    Returns:
        the `sample_group` of the concatenated files.
    """
    f_out = f"{os.path.dirname(samples[0])}/"
    sample_group = re.sub(DATA_DIR, "", f_out).strip("/")
    if not os.path.exists(f"{f_out}merged_1.fq.gz"):
        f_of_sample = " ".join([f"{x}_1.fq.gz" for x in samples])
        cmd1 = f"cat {f_of_sample} > {f_out}merged_1.fq.gz"
        subprocess.check_call(cmd1, shell=True)
        logger.info(f"Concatenated fastq files for {sample_group}_1")

    if not os.path.exists(f"{f_out}merged_2.fq.gz"):
        f_of_sample = " ".join([f"{x}_2.fq.gz" for x in samples])
        cmd2 = f"cat {f_of_sample} > {f_out}merged_2.fq.gz"
        subprocess.check_call(cmd2, shell=True)
        logger.info(f"Concatenated fastq files for {sample_group}_2")
    return sample_group


if __name__ == "__main__":
    logger.info("\x1b[31;1m" + "/*** GDC RNA-Seq pipeline started! ***/" + "\x1b[0m")
    ###################################################################
    #                        Get design matrix                        #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 1: get design matrix" + "\x1b[0m")
    if os.path.exists(f"{DATA_DIR}/design_matrix.csv"):
        design_mat = pd.read_csv(f"{DATA_DIR}/design_matrix.csv")
        logger.info("Read design matrix from file")
    else:
        # Get sample paths for original and new files
        original_samples = ListFastqFiles(group="original")
        new_samples = ListFastqFiles(group="new")
        all_samples = original_samples + new_samples

        # Make design matrix of raw data files
        design_mat = Samples2DesignMatrix(all_samples)
        design_mat.to_csv(f"{DATA_DIR}/design_matrix.csv", index=False)
        logger.info("Created design matrix")

    ###################################################################
    #                 Concatenate multi-lane samples                  #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 2: concatenate multi-lane samples" + "\x1b[0m")
    sm_multi_lanes = design_mat[design_mat["filename"].str.contains(",")]
    filenames = sm_multi_lanes.apply(
        lambda x: [
            f"{DATA_DIR}/{x['sample_group']}/{ele}" for ele in x["filename"].split(",")
        ],
        axis=1,
    ).tolist()
    for f in filenames:
        sg = ConcatSamples(f)
        design_mat.loc[design_mat["sample_group"] == sg, "filename"] = "merged"

    assert not any(design_mat.filename.str.contains(","))

    ###################################################################
    #                     Run fastp on all fq files                   #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 3: QC and preprocess with fastp" + "\x1b[0m")
    filenames = (
        DATA_DIR + "/" + design_mat["sample_group"] + "/" + design_mat["filename"]
    )
    filenames = filenames.tolist()
    sample_groups = [x.replace("/", "_") for x in design_mat["sample_group"]]
    for i, (f, sg) in enumerate(zip(filenames, sample_groups)):
        if os.path.exists(f"{RES_DIR}/fastp/{sg}_fastp.json"):
            continue
        logger.info(f"Running fastp on sample {sg}...")
        subprocess.check_call(
            f"""{FASTP_PATH} -V -i {f}_1.fq.gz -I {f}_2.fq.gz \
                -o {RES_DIR}/data/{sg}_1.fq.gz \
                -O {RES_DIR}/data/{sg}_2.fq.gz \
                --html {RES_DIR}/fastp/{sg}_fastp.html \
                --json {RES_DIR}/fastp/{sg}_fastp.json \
                -w 4""",
            shell=True,
        )
        logger.info(f"Generated fastp report for sample {sg}")

    subprocess.check_call(
        f"{MULTIQC_PATH} {RES_DIR}/fastp/ -m fastp -o {RES_DIR}/fastp/multiqc/",
        shell=True,
    )

    if not os.path.exists(f"{RES_DIR}/fastqc/multiqc_report.html"):
        subprocess.check_call(
            f"{FASTQC_PATH} {RES_DIR}/data/* --noextract -o {RES_DIR}/fastqc/ -t 4"
        )
        subprocess.check_call(
            f"{MULTIQC_PATH} {RES_DIR}/fastp/ -m fastqc -o {RES_DIR}/fastqc/",
            shell=True,
        )

    ###################################################################
    #                 Align sequences and call counts                 #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 4: STAR alignment" + "\x1b[0m")

    # Build the STAR index if it's not already built
    if not os.path.exists(f"{STAR_INDEX_DIR}/Genome"):
        logger.info("STAR index not found. Building now...")
        subprocess.check_call(
            f"""STAR \
                --runMode genomeGenerate \
                --genomeDir {STAR_INDEX_DIR} \
                --genomeFastaFiles {REFERENCE_GENOME_PATH} \
                --sjdbOverhang 100 \
                --sjdbGTFfile {GENCODE_PATH} \
                --runThreadN 8 \
                --outFileNamePrefix {WORK_DIR}/logs/star_index""",
            shell=True,
        )
        logger.info(f"STAR index built to {STAR_INDEX_DIR}")

    # Run STAR for each sample if output files are not found
    sg_sms = design_mat["sample_group"].apply(os.path.basename)
    for i, (sg, sm) in enumerate(zip(sample_groups, sg_sms)):
        if os.path.exists(f"{RES_DIR}/counts/{sg}.tsv") and os.path.exists(
            f"{RES_DIR}/bam/{sg}"
        ):
            continue
        os.makedirs(f"{RES_DIR}/bam/{sg}", exist_ok=True)
        logger.info(f"Aligning sample {sg}")
        subprocess.check_call(
            f"""{STAR_PATH} \
                --readFilesIn {RES_DIR}/data/{sg}_1.fq.gz {RES_DIR}/data/{sg}_2.fq.gz \
                --outSAMattrRGline ID:{sg} SM:{sm} \
                --alignIntronMax 1000000 \
                --alignIntronMin 20 \
                --alignMatesGapMax 1000000 \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 8 \
                --alignSoftClipAtReferenceEnds Yes \
                --chimJunctionOverhangMin 15 \
                --chimMainSegmentMultNmax 1 \
                --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
                --chimSegmentMin 15 \
                --genomeDir {STAR_INDEX_DIR} \
                --genomeLoad NoSharedMemory \
                --limitSjdbInsertNsj 1200000 \
                --outFileNamePrefix {RES_DIR}/bam/{sg}/{sg} \
                --outFilterIntronMotifs None \
                --outFilterMatchNminOverLread 0.33 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverLmax 0.1 \
                --outFilterMultimapNmax 20 \
                --outFilterScoreMinOverLread 0.33 \
                --outFilterType BySJout \
                --outSAMattributes NH HI AS nM NM ch \
                --outSAMstrandField intronMotif \
                --outSAMtype BAM Unsorted \
                --outSAMunmapped Within \
                --quantMode TranscriptomeSAM GeneCounts \
                --readFilesCommand zcat \
                --runThreadN 8 \
                --twopassMode Basic""",
            shell=True,
        )
        shutil.move(
            f"{RES_DIR}/bam/{sg}/{sg}ReadsPerGene.out.tab", f"{RES_DIR}/counts/{sg}.tsv"
        )
        logger.info(f"Counts for sample {sg} generated")

    ###################################################################
    #                      Get TPM using Salmon                       #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 5: Get TPM values using Salmon" + "\x1b[0m")
    for sg in sample_groups:
        if os.path.exists(f"{RES_DIR}/tpm/{sg}"):
            continue
        logger.info(f"Calling TPM values for sample {sg}")
        subprocess.check_call(
            f"""{SALMON_PATH} quant \
                -t {REFERENCE_TRANSCRIPTOME_PATH} \
                -l A \
                -a {RES_DIR}/bam/{sg}/{sg}Aligned.toTranscriptome.out.bam \
                -o {RES_DIR}/tpm/{sg} \
                --gencode""",
            shell=True,
        )

    ###################################################################
    #                 Combine counts and TPM tables                   #
    ###################################################################
    logger.info("\x1b[33;21m" + "Step 6: Combine counts and TPM tables" + "\x1b[0m")
    counts_table = []
    tpm_table = []
    for sg in sample_groups:
        # Combine counts files into one table
        df = pd.read_table(f"{RES_DIR}/counts/{sg}.tsv", header=None).iloc[:, 0:2]
        df.columns = ["Ensembl", "count"]
        df["sample_group"] = sg
        counts_table.append(df)

        # Combine TPM files, also copy the quant.sf files to a separate folder
        df = pd.read_table(f"{RES_DIR}/tpm/{sg}/quant.sf")[["Name", "TPM"]]
        df.columns = ["Ensembl", "TPM"]
        df["sample_group"] = sg
        tpm_table.append(df)

    counts_table = pd.concat(counts_table, axis=0, ignore_index=True)
    counts_table = counts_table.pivot(
        index="Ensembl", columns="sample_group", values="count"
    )
    tpm_table = pd.concat(tpm_table, axis=0, ignore_index=True)
    tpm_table = tpm_table.pivot(index="Ensembl", columns="sample_group", values="TPM")

    counts_table.to_csv(f"{RES_DIR}/counts.csv")
    tpm_table.to_csv(f"{RES_DIR}/TPM.csv")

    logger.info("\x1b[31;1m" + "Transcript quantification finished!" + "\x1b[0m")
    # Cleanup
    # shutil.rmtree(DATA_DIR)
