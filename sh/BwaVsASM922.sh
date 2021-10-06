#!/usr/bin/env bash
#SBATCH -J bwa-mem
#SBATCH --account=project_2001175
#SBATCH --time 0-00:15:00
#SBATCH -p test
#SBATCH --mem 40GB
#SBATCH -n 20

#samtools faidx Pf_ANC_contigs.fa
#bwa index Pf_ANC_contigs.fa

SAMPLE=${1}

CONTAMS="/projappl/project_2001175/software/bbmap/v38.75/resources/adapters.fa"
bbduk.sh ow=t ftm=5 ref=${CONTAMS} ktrim=r k=23 mink=11 hdist=1 hdist2=1 tbo tpe in=../rawreads/${SAMPLE}_R1.fastq.gz in2=../rawreads/${SAMPLE}_R2.fastq.gz out=stdout.fastq | bbduk.sh int=t ref=phix k=31 hdist=1 in=stdin.fastq out=stdout.fastq | bbduk.sh ow=t int=t qtrim=r trimq=10 in=stdin.fastq out=../qcreads/${SAMPLE}_R1.fastq.gz out2=../qcreads/${SAMPLE}_R2.fastq.gz

bwa mem \
-t 16 \
-R "@RG\tID:0\tSM:Pf_ANC\tLB:FIMM\tPU:Illumina" \
../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.fasta \
../qcreads/${SAMPLE}_R1.fastq.gz ../qcreads/${SAMPLE}_R2.fastq.gz | \
samtools view -bh | \
samtools sort -@ 4 -o ../mapped/${SAMPLE}.ASM922.bwa-mem.bam -

cd ../mapped

samtools index ${SAMPLE}.ASM922.bwa-mem.bam

