#!/usr/bin/env bash
#SBATCH -J shovill
#SBATCH --account=project_2001175
#SBATCH --time 0-07:15:00
#SBATCH -p small
#SBATCH --mem 40GB
#SBATCH -n 20
#SBATCH --gres=nvme:100

samtools faidx Pf_ANC_contigs.fa
bwa index Pf_ANC_contigs.fa

SAMPLE=${1}

CONTAMS="/projappl/project_2001175/software/bbmap/v38.75/resources/adapters.fa"
bbduk.sh ow=t ftm=5 ref=${CONTAMS} ktrim=r k=23 mink=11 hdist=1 hdist2=1 tbo tpe in=../rawreads/${SAMPLE}_R1.fastq.gz in2=../rawreads/${SAMPLE}_R2.fastq.gz out=stdout.fastq | bbduk.sh int=t ref=phix k=31 hdist=1 in=stdin.fastq out=stdout.fastq | bbduk.sh ow=t int=t qtrim=r trimq=10 in=stdin.fastq out=${LOCAL_SCRATCH}/${SAMPLE}_R1.fastq.gz out2=${LOCAL_SCRATCH}/${SAMPLE}_R2.fastq.gz

bwa mem \
-t 16 \
-R "@RG\tID:0\tSM:Pf_ANC\tLB:FIMM\tPU:Illumina" \
../Pf_ANC_assembly/Pf_ANC_contigs.fa \
${LOCAL_SCRATCH}/${SAMPLE}_R1.fastq.gz ${LOCAL_SCRATCH}/${SAMPLE}_R2.fastq.gz | \
samtools view -bh | \
samtools sort -@ 4 -o ../mapped/${SAMPLE}.Pf_ANC.bwa-mem.bam -

mv ${LOCAL_SCRATCH}/${SAMPLE}_R1.fastq.gz ../qcreads
mv ${LOCAL_SCRATCH}/${SAMPLE}_R2.fastq.gz ../qcreads

cd ../mapped

samtools index ${SAMPLE}.Pf_ANC.bwa-mem.bam

