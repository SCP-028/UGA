#!/bin/sh
ROOTDIR=$(pwd)
SOFTWAREDIR="$ROOTDIR/BINF8441/software"
mkdir -p "$SOFTWAREDIR"
cd "$SOFTWAREDIR"

# sratoolkit: fastq file download
if [ ! -d sratoolkit* ]; then
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-ubuntu64.tar.gz
    tar -xzvf sratoolkit.2.9.0-ubuntu64.tar.gz && rm sratoolkit.2.9.0-ubuntu64.tar.gz
fi
export PATH="$SOFTWAREDIR/sratoolkit.2.9.0-ubuntu64/bin:$PATH"

# FastQC: pre-alignment QC
if [ ! -d FastQC ]; then
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip fastqc_v0.11.7.zip && rm fastqc_v0.11.7.zip
fi
export PATH="$SOFTWAREDIR/FastQC:$PATH"

# skewer: trim adapter sequences
# wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz
# tar -xzvf 0.2.2.tar.gz && rm 0.2.2.tar.gz
# cd skewer-0.2.2
# mkdir bin
# sed -i "s|/usr/local/bin|${ROOTDIR}/BINF8441/software/skewer-0.2.2/bin|" Makefile
# make
# make install
# export PATH="$SOFTWAREDIR/skewer-0.2.2/bin:$PATH"
# cd "$SOFTWAREDIR"

# STAR: read aligner
if [! -d STAR-2.6.0a ]; then
    wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
    tar -xzvf 2.6.0a.tar.gz && rm 2.6.0a.tar.gz
fi
export PATH="$SOFTWAREDIR/STAR-2.6.0a/bin/Linux_x86_64:$PATH"

# Samtools: view BAM files
if [ ! -d samtools-1.8/pkg ]; then
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    bzip2 -d samtools-1.8.tar.bz2
    tar -xvf samtools-1.8.tar && rm samtools-1.8.tar
    cd samtools-1.8
    mkdir pkg
    autoheader
    autoconf -Wno-syntax
    ./configure --prefix="$SOFTWAREDIR/samtools-1.8/pkg"
    make -j5
    make install
    cd "$SOFTWAREDIR"
fi
export PATH="$SOFTWAREDIR/samtools-1.8/pkg/bin:$PATH"

# HTSeq: BAM file to counts
pip install HTSeq

# subread: BAM file to counts
if [ ! -d subread* ]; then
    wget https://downloads.sourceforge.net/project/subread/subread-1.6.1/subread-1.6.1-Linux-x86_64.tar.gz
    tar -xzvf subread-1.6.1-Linux-x86_64.tar.gz && rm subread-1.6.1-Linux-x86_64.tar.gz
fi
export PATH="$SOFTWAREDIR/subread-1.6.1-Linux-x86_64/bin:$PATH"

# download fastq files
cd "$ROOTDIR/BINF8441"
FASTQDIR="$ROOTDIR/BINF8441/fastq"
mkdir "$FASTQDIR"
cd "$FASTQDIR"
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2650363
# --split-3 because these used paired-end sequencing
prefetch SRX2884433
fastq-dump --split-3 "$HOME/ncbi/public/sra/SRX2884433.sra"  # control iron 1

# pre-alignment quality control reports
mkdir qc
fastqc SRX2884433_1.fastq SRX2884433_2.fastq -o "./qc" -t 2
# Check quality control reports before moving forward! #

################################################################################
# The following step is skipped because the produced files can't be fed to STAR,
# and the QC report didn't find any overrepresented sequences.
# remove adapter sequences
# skewer \
#     -n \
#     -m pe \
#     -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
#     SRX2884433_1.fastq SRX2884433_2.fastq
################################################################################

mkdir log raw_sequence
mv *.log ./log/
# mkdir trimmed_sequence
# mv *trimmed.* ./trimmed_sequence/
mv *.fastq ./raw_sequence/

# generate genome index files
if [ ! -d annotation ]; then
    mkdir annotation
    cd annotation
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    gunzip gencode.v28.annotation.gtf.gz
    cd "$FASTQDIR"
fi

if [ ! -d "$FASTQDIR/GenomeDir"]; then
    mkdir -p GenomeDir
    STAR \
        --runThreadN 30 \
        --runMode genomeGenerate \
        --genomeDir "$FASTQDIR/GenomeDir" \
        --genomeFastaFiles "$FASTQDIR/annotation/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" \
        --sjdbGTFfile "$FASTQDIR/annotation/gencode.v28.annotation.gtf" \
        --sjdbOverhang 100
    mv Log.out "$FASTQDIR/log/genome.log.out"
fi

# map reads to the genome
mkdir -p ctrlGenome aligned
## alignment 1st pass
STAR \
    --runThreadN 30 \
    --genomeDir "$FASTQDIR/GenomeDir" \
    --readFilesIn "$FASTQDIR/raw_sequence/SRX2884433_1.fastq" "$FASTQDIR/raw_sequence/SRX2884433_2.fastq" \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMtype None \
    --outSAMmode None
mv SJ.out.tab "$FASTQDIR/annotation/ctrl.SJ.out.tab"
## intermediate index generation
STAR \
    --runThreadN 30 \
    --runMode genomeGenerate \
    --genomeDir "$FASTQDIR/ctrlGenome" \
    --genomeFastaFiles GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --sjdbOverhang 100 \
    --sjdbFileChrStartEnd "$FASTQDIR/annotation/ctrl.SJ.out.tab"
## alignment 2nd pass
STAR \
    --runThreadN 20 \
    --genomeDir "$FASTQDIR/ctrlGenome" \
    --readFilesIn "$FASTQDIR/raw_sequence/SRX2884433_1.fastq" "$FASTQDIR/raw_sequence/SRX2884433_2.fastq" \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 0 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMheaderHD @HD VN:1.4
mv Log.out "$FASTQDIR/log/ctrl.align.log.out"
mv Log.progress.out "$FASTQDIR/log/ctrl.align.log.progress.out"
mv Log.final.out "$FASTQDIR/log/ctrl.align.log.final.out"
mv Aligned.sortedByCoord.out.bam "$FASTQDIR/aligned/ctrl.aligned.bam"

# convert BAM to SAM
cd aligned
samtools view -F 4 -@ 20 -o ctrl.aligned.sam ctrl.aligned.bam

# aligned SAM file to gene count (TCGA pipeline)
htseq-count \
-m intersection-nonempty \
-i gene_id \
-r pos \
-s no \
ctrl.aligned.sam "$FASTQDIR/annotation/gencode.v28.annotation.gtf" |
tee ctrl.counts.htseq

## alternatively, use featureCounts (should be much faster)
featureCounts \
-p \
-g gene_id \
-T 20 \
-a "$FASTQDIR/annotation/gencode.v28.annotation.gtf" \
-o ctrl.counts.featureCounts \
ctrl.aligned.sam
