#!/bin/bash
PROJECT_NAME="GSE75168" # single-ended
declare -a SRX=("SRX1438068" "SRX1438069" "SRX1438070" "SRX1438074" "SRX1438075" "SRX1438076")
# MCF10A - 1, 2, 3    |  MDA MB-231 - 1, 2, 3
ROOTDIR=$(pwd)
SOFTWAREDIR="$ROOTDIR/software"
FASTQDIR="$ROOTDIR/fastq"
USE_CACHE=true  # if false, remove all existing data / software

###############################################################################
#                Shouldn't need to modify below this line                     #
###############################################################################
mkdir -p "$FASTQDIR/qc" "$SOFTWAREDIR" && cd "$SOFTWAREDIR"

# sratoolkit: fastq file download
if [ ! -d sratoolkit* ] || [ ! USE_CACHE ]; then
    rm -rf ./sratoolkit*
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -xzvf sratoolkit.current-ubuntu64.tar.gz && rm sratoolkit.current-ubuntu64.tar.gz
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'sratoolkit*')/bin:$PATH"

# FastQC: pre-alignment QC
if [ ! -d FastQC* ]|| [ ! USE_CACHE ]; then
    rm -rf ./FastQC
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip fastqc_v0.11.7.zip && rm fastqc_v0.11.7.zip
    chmod +x FastQC/fastqc
fi
export PATH="$SOFTWAREDIR/FastQC:$PATH"

if [ ! -d Trimmomatic* ] || [ ! USE_CACHE]; then
    wget "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip"
    unzip Trimmomatic-0.38.zip && rm Trimmomatic-0.38.zip
fi

# HISAT2: alignment
if [ ! -d hisat2* ]|| [ ! USE_CACHE ]; then
    rm -rf ./hisat2*
    wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip"
    unzip hisat2-2.1.0-Linux_x86_64.zip && rm hisat2-2.1.0-Linux_x86_64.zip
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'hisat2*' | tail -n 1):$PATH"

# Samtools: view BAM files
if [ ! -d samtools* ]|| [ ! USE_CACHE ]; then
    rm -rf ./samtools*
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
if [ ! -d stringtie* ]|| [ ! USE_CACHE ]; then
    rm -rf ./stringtie*
    wget "http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz"
    tar -xzvf stringtie-1.3.4d.Linux_x86_64.tar.gz && rm stringtie-1.3.4d.Linux_x86_64.tar.gz
fi
export PATH="$(find $SOFTWAREDIR -type d -name 'stringtie*'):$PATH"

# download genome index files
cd $ROOTDIR
if [ ! -d "$ROOTDIR/grch38" ]|| [ ! USE_CACHE ];then
    rm -rf "$ROOTDIR/grch38"
    wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz"
    tar xzvf grch38.tar.gz && rm grch38.tar.gz
    bash grch38/make_grch38.sh
    wget "ftp://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.gff3.gz"
    gunzip Homo_sapiens.GRCh38.94.gff3.gz
    mv Homo_sapiens.GRCh38.94.gff3 ./grch38/
fi

#####################################################

mkdir -p "$FASTQDIR/$PROJECT_NAME" && cd "$FASTQDIR/$PROJECT_NAME"
if [ ! USE_CACHE ];then
    rm -rf $HOME/ncbi/public/sra
fi
for f in "$SRX[@]"
do
    # download fastq files
    if [ ! -f "$f.fastq" ];then
        prefetch "$f" && vdb-dump -f fastq "$HOME/ncbi/public/sra/$f.sra" > "$f.fastq"
    fi
done

# pre-alignment quality control reports
if [ ! -d "$FASTQDIR/$PROJECT_NAME/qc/before_trim" ];then
    mkdir -p "$FASTQDIR/$PROJECT_NAME/qc/before_trim"
    fastqc $(printf "%s.fastq " "${SRX[@]}") -o "$FASTQDIR/$PROJECT_NAME/qc/before_trim" -t ${#SRX[@]}
fi

for f in "$SRX[@]"
do
    # trim adapters
    if [ ! -f "$FASTQDIR/$PROJECT_NAME/$f.trimmed.fq" ];then
        java -jar "$SOFTWAREDIR/Trimmomatic-0.38/trimmomatic-0.38.jar" SE -threads 30 -phred33 "$f.fastq" "$f.trimmed.fq" \
            ILLUMINACLIP:"$SOFTWAREDIR/Trimmomatic-0.38/adapters/TruSeq3-SE.fa":2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && rm "$f.fastq"
    fi
done

if [ ! -d "$FASTQDIR/$PROJECT_NAME/qc/after_trim" ];then
    mkdir -p "$FASTQDIR/$PROJECT_NAME/qc/after_trim"
    fastqc $(printf "%s.trimmed.fq " "${SRX[@]}") -o "$FASTQDIR/$PROJECT_NAME/qc/after_trim/" -t ${#SRX[@]}
fi

# Check quality control reports before moving forward!
###########################################################################
# Per base sequence quality should be mostly in the green range           #
# Per sequence GC content should be within the 30%-70% range              #
# Sequence length distribution could help determine MINLEN in trimmomatic #
# Adapter content should be around 0 after trimming                       #
###########################################################################

for f in "$SRX[@]"
do
    # align reads
    hisat2 --dta -p 30 \
        -x "$ROOTDIR/grch38/genome" \
        -U "$FASTQDIR/$PROJECT_NAME/$f.trimmed.fq" \
        -S "$f.sam" && rm "$f.trimmed.fq"
    # convert SAM to BAM
    samtools sort -@ 30 -o "$f.bam" "$f.sam" && rm "$f.sam"
    # transcript assembly
    stringtie -p 30 -G "$ROOTDIR/grch38/Homo_sapiens.GRCh38.94.gff3" -o "$f.gtf" -l "$f" "$f.bam"
done

ls | grep "gtf$" > mergelist.txt
stringtie --merge -p 30 -G "$ROOTDIR/grch38/Homo_sapiens.GRCh38.94.gff3" -o stringtie_merged.gtf mergelist.txt

for f in "$SRX[@]"
do
    rm "$f.gtf"
    mkdir -p "$FASTQDIR/$PROJECT_NAME/ballgown/$f"
    stringtie -e -B -p 30 -G stringtie_merged.gtf -o "$FASTQDIR/$PROJECT_NAME/ballgown/$f/$f.gtf" "$f.bam"
done
cd $ROOTDIR

# Use R package ballgown for differential expression analysis
