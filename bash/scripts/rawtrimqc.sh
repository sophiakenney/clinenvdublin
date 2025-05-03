#!/bin/bash
# Raw Quality Control

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate bioinfo

# ---- Run FastQC and MultiQC on Raw Reads ----

cd /your/path/to/rawreads/0.0-raw

# --- run fastqc
fastqc *.fq.gz -o "/your/path/to/rawqc_out/0.1-rawqc" -t 10

# --- multiqc

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

multiqc /your/path/to/rawqc_out/0.1-rawqc/*_fastqc.zip --interactive

# --- run seqkit
seqkit stats *.gz -T > /your/path/to/rawqc_out/0.1-rawqc/rawstats.txt


# ---- Trim Raw Reads with Trimmomatic ----

# --- set paths for trimmomatic and trim qc
RAWDIR=/your/path/to/rawreads/0.0-raw
OUTDIR=/your/path/to/trimmed_out/1.0-trimmed
ADAPTERS=/your/path/to/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/

cd ${RAWDIR}

for f in *_1.fq.gz
do
outfile=${f%_1.fq.gz}.fq.gz
echo "${f}"

/your/path/to/trimmomatic PE "${f}" "${f%_1.fq.gz}_2.fq.gz" \
-baseout "${OUTDIR}/${outfile}" -threads 10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ILLUMINACLIP:${ADAPTERS}TruSeq3-PE-2.fa:2:30:10

done

# ---- Trimmed Read QC ----

cd ${OUTDIR}

# --- run fastqc
fastqc *P.fq.gz -o "/your/path/to/trimqc_out/1.1-trimqc" -t 10

# --- multiqc

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

multiqc /your/path/to/trimqc_out/1.1-trimqc/*_fastqc.zip \
    -o  /your/path/to/trimqc_out/1.1-trimqc \
    --interactive

# --- run seqkit
seqkit stats *.gz -T > /your/path/to/trimqc_out/trimstats.txt

# --- write job end time
echo "Job ended at $(date)"
