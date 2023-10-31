#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=cull0084@umn.edu 
#SBATCH --job-name jg_snake
#SBATCH -o ./Jobs/%j.pile_gsnake.out
#SBATCH -e ./Jobs/%j.pile_gsnake.err
#SBATCH -p amdsmall,amdlarge,amd512,small,large,ram256g

set -e

conda activate snake718
cd $SLURM_SUBMIT_DIR

BAM_DIR=/scratch.global/friedlab_LOW_PILEUP
PHASED_POP=/panfs/jay/groups/0/fried255/fried255/working/pipeline
SUB_POP=/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass

snakemake -s slopi_snake.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$BAM_DIR,$PHASED_POP,$SUB_POP" \
    --profile slurm.go_wags/ \
    --configfile config.yaml \
    --rerun-triggers mtime \
    --keep-going

