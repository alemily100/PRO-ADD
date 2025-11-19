#!/bin/bash
#SBATCH --job-name=array-mayo
#SBATCH --output=mayo-%a.txt
#SBATCH --error=mayo-%a.txt
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-8
#SBATCH --mail-user=emily.alger@icr.ac.uk
#SBATCH --mail-type=ALL

source ~/.bashrc
eval "$(mamba shell hook --shell bash)"
mamba activate mamba-emily

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun Rscript /home/ealger/pro_add/sim_array_hypothetical_partial.R $SLURM_ARRAY_TASK_ID
