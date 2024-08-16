import os
import itertools
import pandas as pd
from pathlib import Path

singularity: config['sif']

rule all:
    input:
       ## begin main genotyping and imputation
        expand(
            "results/{ref}/impute/{breed}/combine/{breed}.dr2_fltr_imputed.snps.{ref}.{date}.vcf.gz.covtsf",
            ref=config['refs'],
            breed=config['breeds'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M OR Y
            date=config['date'],
        )

include: "rules/genotype.smk"
include: "rules/impute.smk"

