import os
import itertools
import pandas as pd
from pathlib import Path

localrules: gather_gatk_summaries,
            gather_imputed_gatk_summaries

singularity: config['sif']
include: "src/utils.py"

REFS = ["cf4"]

rule all:
    input:
       ## begin main genotyping and imputation
        expand(
            "results/{ref}/vcf/{breed}/combine/{breed}.fltr_cyrano_imputed.snps.{ref}.20231018.vcf.{ext}",
            ref=REFS,
            breed=config["breeds"],
           #chrom=['chr25'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M CORRECT?
            ext=["gz","gz.tbi"]
        ),
       #expand(
       #    "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
       #    ref=REFS,
       #    wga=WGA,
       #    nowga=NOWGA
       #   #ext=["gz","gz.tbi"]
       #),
       #"results/imputed_subset_concord.gatk.csv"
       # BESPOKE
        expand(
             "results/{ref}/vcf/gtdn/filter/gtdn.fltr_gtdn_imputed_chr25.snps.{ref}.20231018.vcf.gz",
             ref=REFS,
         ),
        expand(
             "results/{ref}/vcf/gshp/filter/gshp.fltr_gshp_imputed_chr24.snps.{ref}.20231018.vcf.gz",
             ref=REFS,
         ),


            
            

#include: "rules/subset_concord.smk"
include: "rules/genotype.smk"
include: "rules/impute.smk"
#include: "rules/impute_concord.smk"
include: "rules/bespoke_impute.smk"

