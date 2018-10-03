#!/bin/sh
# GSE78011
declare -a SRX=("SRX1589783" "SRX1589784" "SRX1589785" "SRX1589789" "SRX1589790")
# MCF7 - 1, 2, 3    |  MDA MB-231 - 1, 2
ROOTDIR=$(pwd)
SOFTWAREDIR="$ROOTDIR/software"
FASTQDIR="$ROOTDIR/fastq"
mkdir -p "$FASTQDIR/qc" "$SOFTWAREDIR"

cd "$SOFTWAREDIR"

# sratoolkit: fastq file download
if [ ! -d sratoolkit* ]; then
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzvf sratoolkit.current-ubuntu64.tar.gz && rm sratoolkit.current-ubuntu64.tar.gz
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'sratoolkit*')/bin:$PATH"

# FastQC: pre-alignment QC
if [ ! -d FastQC* ]; then
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip fastqc_v0.11.7.zip && rm fastqc_v0.11.7.zip
fi
export PATH="$SOFTWAREDIR/FastQC:$PATH"

# HISAT2: alignment
if [ ! -d hisat2* ]; then
    wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip"
    unzip hisat2-2.1.0-Linux_x86_64.zip && rm hisat2-2.1.0-Linux_x86_64.zip
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'hisat2*' | tail -n 1):$PATH"

# Samtools: view BAM files
if [ ! -d samtools* ]; then
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    bzip2 -d samtools-1.9.tar.bz2 && tar -xvf samtools-1.9.tar && rm samtools-1.9.tar
    cd samtools-1.9
    mkdir pkg
    autoheader
    autoconf -Wno-syntax
    ./configure --prefix="$SOFTWAREDIR/samtools-1.9/pkg"
    make -j5
    make install
    cd "$SOFTWAREDIR"
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'samtools*')/pkg/bin:$PATH"

# StringTie: transcript assembler
if [ ! -d stringtie* ]; then
    wget "http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz"
    tar -xzvf stringtie-1.3.4d.Linux_x86_64.tar.gz && rm stringtie-1.3.4d.Linux_x86_64.tar.gz
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'stringtie*'):$PATH"

# download genome index files
if [ ! -d "$ROOTDIR/grch38" ];then
    wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz"
    tar xzvf grch38.tar.gz && rm grch38.tar.gz
    bash grch38/make_grch38.sh
    wget "ftp://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.gff3.gz"
    gunzip Homo_sapiens.GRCh38.94.gff3.gz
    mv Homo_sapiens.GRCh38.94.gff3 ./grch38/
fi

#####################################################

# download fastq files
cd "$FASTQDIR"
for f in "$SRX[@]"
do
    prefetch "$f"
    vdb-dump -f fastq "$HOME/ncbi/public/sra/$f.sra" > "$f.fastq"
    # pre-alignment quality control reports
    # mkdir -p "$FASTQDIR/qc/$f"
    # fastqc "$f.fastq" -o "./qc/$f" -t 4
    ########################################################
    # Check quality control reports before moving forward! #
    ########################################################
    # align reads
    hisat2 --dta -p 30 -x "$ROOTDIR/grch38/genome" -U "$FASTQDIR/$f.fastq" -S "$f.sam" && rm "$f.fastq"
    # convert SAM to BAM
    # samtools view -@ 20 -b -o "$f.sam" "$f.bam"
    samtools sort -@ 30 -o "$f.bam" "$f.sam" && rm "$f.sam"
    # transcript assembly
    stringtie -p 30 -G "$ROOTDIR/grch38/Homo_sapiens.GRCh38.94.gff3" -o "$f.gtf" -l "$f" "$f.bam"
done

ls | grep "gtf$" > mergelist.txt
stringtie --merge -p 30 -G "$ROOTDIR/grch38/Homo_sapiens.GRCh38.94.gff3" -o stringtie_merged.gtf mergelist.txt

for f in "$SRX[@]"
do
    mkdir -p "$FASTQDIR/ballgown/$f"
    stringtie -e -B -p 30 -G stringtie_merged.gtf -o "$FASTQDIR/ballgown/$f/$f.gtf" "$f.bam"
done

# Use R package ballgown for differential expression analysis
