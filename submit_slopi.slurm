#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=cull0084@umn.edu 
#SBATCH --job-name jg_snake
#SBATCH -o ./Jobs/%j.pile_gsnake.out
#SBATCH -e ./Jobs/%j.pile_gsnake.err
#SBATCH -p msismall,msilarge

set -e

conda activate snake718
cd $SLURM_SUBMIT_DIR

BAM3_DIR=/scratch.global/friedlab_NEO_CF3
BAM4_DIR=/scratch.global/friedlab_NEO_CF4
PHASED_POP=/panfs/jay/groups/0/fried255/fried255/working/pipeline

snakemake -s slopi_snake.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$BAM3_DIR,$BAM4_DIR,$PHASED_POP" \
    --profile slurm.go_wags/ \
    --configfile config.yaml \
    --rerun-triggers mtime \
    --keep-going

