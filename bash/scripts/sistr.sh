#!/bin/bash
# Run SISTR

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate sistr

#set directory paths
ASSEMB=/your/path/to/assemblies/3.0-unicycler_entero/assemblies
OUT=/our/path/to/3.2-sero/sistr

#change to working directory
cd ${ASSEMB}

for f in *_assembly.fasta ;
do

  echo ${f%_assembly.fasta}

  sistr --qc \
  -vv \
  -t 4 \
  --alleles-output ${f%_assembly.fasta}_allele-results.json \
  --novel-alleles ${f%_assembly.fasta}_novel-alleles.fasta \
  --cgmlst-profiles ${f%_assembly.fasta}_cgmlst-profiles.csv \
  -f tab \
  -o ${f%_assembly.fasta}-output.tab ${f} ;

done

# move results to output directory:

# create subdirectories if they don't already exist

mkdir -p ${OUT}/alleles
mkdir -p ${OUT}/novel
mkdir -p ${OUT}/cgmlst
mkdir -p ${OUT}/output

# move files

for f in *_assembly.fasta

do
        mv ${f%_assembly.fasta}_allele-results.json ${OUT}/alleles
        mv ${f%_assembly.fasta}_novel-alleles.fasta ${OUT}/novel
        mv ${f%_assembly.fasta}_cgmlst-profiles.csv ${OUT}/cgmlst
        mv ${f%_assembly.fasta}-output.tab ${OUT}/output

done


#create one file

cd ${OUT}/output

{ head -n 1 ARS1-output.tab; tail -n +2 -q *.tab; } > sistr_output.txt

cd ../cgmlst

{ head -n 1 ARS1_cgmlst-profiles.csv; tail -n +2 -q *profiles.csv; } > sistr_cgmlst.csv

# --- write job end time
echo "Job ended at $(date)"
