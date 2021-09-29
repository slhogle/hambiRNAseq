#!/usr/bin/env bash
#SBATCH -J shovill
#SBATCH --account=project_2001175
#SBATCH --time 0-07:15:00
#SBATCH -p small
#SBATCH --mem 40GB
#SBATCH -n 8
#SBATCH --gres=nvme:100

# needed to make singularity work
module purge

export SING_IMAGE=/projappl/project_2001175/singularity/torstyverse.sif

singularity_wrapper exec shovill \
--R1 ../rawreads/Pf_ANC_R1.fastq.gz \
--R2 ../rawreads/Pf_ANC_R2.fastq.gz \
--outdir ../Pf_ANC_assembly  \
--gsize 6.7M \
--cpus 8 \
--ram 40 \
--trim \
--tmpdir ${LOCAL_SCRATCH}
