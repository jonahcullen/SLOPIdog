#!/bin/bash -l
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=
#SBATCH --job-name jg_snake
#SBATCH -o ./Jobs/%j.pile_gsnake.out
#SBATCH -e ./Jobs/%j.pile_gsnake.err
#SBATCH -p 

set -e
source ~/.bashrc

conda activate snakemake
cd $SLURM_SUBMIT_DIR

BAM_DIR=/path/to/bam/directory
PHASED_POP=/path/to/reference_panel/directory
SUB_POP=/path/to/auxillary/directory/

snakemake -s slopi_snake.smk \
    --use-singularity \
    --singularity-args "-B $PWD,$BAM_DIR,$PHASED_POP,$SUB_POP" \
    --profile slurm.go_wags/ \
    --configfile config.yaml \
    --rerun-triggers mtime \
    --keep-going

