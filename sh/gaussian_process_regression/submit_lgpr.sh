#!/bin/bash -l
#SBATCH --job-name=lgpr
#SBATCH --account=project_2005776
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=10:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000

# Load r-env-singularity
module load r-env-singularity/4.0.2

# Clean up prior TMPDIR and OMP_NUM_THREADS from .Renviron
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2005776" >> ~/.Renviron

# Match thread and core numbers
echo "OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save $1
