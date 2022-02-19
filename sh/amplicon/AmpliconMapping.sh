#!/usr/bin/env bash
#SBATCH -J iina_map
#SBATCH --time 0-15:30:00
#SBATCH -p serial
#SBATCH --mem 32GB
#SBATCH -n 8

REF="HAMBI-30.16S-amplicon-mapdb.ffn"

cat $1 | while read MYLIB;
do

mkdir fastqc_map/${MYLIB}

# map to HAMBI 16S sequences
bbmap.sh ow=t int=f maxindel=20 minid=0.85 ambiguous=best \
ref=${REF} \
ehist=fastqc_map/${MYLIB}/${MYLIB}.ehist \
qahist=fastqc_map/${MYLIB}/${MYLIB}.qahist \
mhist=fastqc_map/${MYLIB}/${MYLIB}.mhist \
idhist=fastqc_map/${MYLIB}/${MYLIB}.idhist \
scafstats=fastqc_map/${MYLIB}/${MYLIB}.scafstats \
statsfile=fastqc_map/${MYLIB}/${MYLIB}.mapstats \
in=processedreads/${MYLIB}.06-maxee.fastq.gz \
out=mappedreads/${MYLIB}.sam

samtools stats mappedreads/${MYLIB}.sam > fastqc_map/${MYLIB}/${MYLIB}.samtoolsstats

# get number of reads mapped to each 16S sequence
pileup.sh ow=t in=mappedreads/${MYLIB}.sam \
out=mappedreads/${MYLIB}.coverage \
rpkm=mappedreads/${MYLIB}.rpkm

done
