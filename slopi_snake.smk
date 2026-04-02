import os
import itertools
import pandas as pd
from pathlib import Path

singularity: config['sif']

rule all:
    input:
        glimpse_stages = expand("results/glimpse/{ref}/{ref}.{stage}.glimpse_imputed.vcf.gz",
            ref = config['refs'],
            stage = ["unfiltered", "filtered"])

include: "rules/genotype.smk"
include: "rules/glimpse_process.smk"