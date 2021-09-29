#!/usr/bin/env bash
#SBATCH -J octopus-indv
#SBATCH --account=project_2001175
#SBATCH --time 0-10:15:00
#SBATCH -p small
#SBATCH --mem 18GB
#SBATCH -n 8

# --gres=nvme:100

# needed to make singularity work
module purge

export SING_IMAGE=/projappl/project_2001175/singularity/octopus.sif

FOREST=/projappl/project_2001175/software/octopus/germline.v0.7.4.forest

SAMPLE=${1}

singularity_wrapper exec octopus \
-R ../GCF_000009225.2_SBW25/GCF_000009225.2_ASM922v1_genomic.fasta \
-I ../mapped/${SAMPLE}.ASM922.bwa-mem.bam \
-C individual \
--organism-ploidy 1 \
--disable-downsampling \
--bamout ../calls/${SAMPLE}.ASM922.bwa-mem.octopus.bam \
--sequence-error-model PCR.NOVASEQ \
--forest ${FOREST} \
-o ../calls/${SAMPLE}.ASM922.bwa-mem.octopus.vcf.gz \
--threads 8 \
--target-read-buffer-memory=16GB \
--temp-directory-prefix ../calls/${SAMPLE}-temp
