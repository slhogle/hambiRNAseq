#!/usr/bin/env bash
#SBATCH -J mapcount
#SBATCH --time 0-5:00:00
#SBATCH -p serial
#SBATCH --mem 80GB
#SBATCH -n 8

MYLIB=$1

bbmap.sh \
in=processedreads/${MYLIB}.05-mRNA.fastq.gz \
ref=genomes/whole_community_assemblies.fna \
ambiguous=all \
ehist=fastqc_map/${MYLIB}.ehist \
qahist=fastqc_map/${MYLIB}.qahist \
mhist=fastqc_map/${MYLIB}.mhist \
idhist=fastqc_map/${MYLIB}.idhist \
scafstats=fastqc_map/${MYLIB}.scafstats \
statsfile=fastqc_map/${MYLIB}.mapstats \
out=mappedreads/${MYLIB}.sam

samtools stats mappedreads/${MYLIB}.sam > fastqc_map/${MYLIB}.samtoolsstats

cat mappedreads/${MYLIB}.sam | samtools sort -O SAM -o mappedreads/${MYLIB}.sorted.sam

featureCounts -T 2 -F SAF -O -M -p -P -B -g ID \
-a genomes/whole_community_assemblies.gaf \
-o mappedreads/${MYLIB}.feats.tsv \
mappedreads/${MYLIB}.sorted.sam
