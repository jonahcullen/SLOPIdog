import os
import itertools
import pandas as pd
from pathlib import Path

localrules: gather_gatk_summaries,
            imputed_list,
            gather_imputed_gatk_summaries

singularity: config['sif']
include: "src/utils.py"

REFS = ["cf3","cf4"]
#SAMPLES = [
#    "D007436_wga","D007436_nowga",
#    "D007619_wga","D007619_nowga",
#    "D007974_wga","D007974_nowga",
#    "D007986_wga","D007986_nowga"
#]

WGA = ["D007436_wga","D007619_wga","D007974_wga","D007986_wga"]
NOWGA = ["D007436_nowga","D007619_nowga","D007974_nowga","D007986_nowga"]


rule all:
    input:
       ## sample VCF
       #expand(
       #   #"results/{ref}/vcf/subset_lowpass.{ref}.vcf.{ext}",
       #    "results/{ref}/vcf/{sample}/{sample}_lowpass.{ref}.vcf.{ext}",
       #    ref=REFS,
       #    sample=WGA+NOWGA,
       #    ext=["gz","gz.tbi"]
       #),
       ## pairwise bcftools stats
       #expand(
       #    "results/{ref}/concordance/bcftools/{wga}_x_{nowga}_lowpass.{ref}.stats.txt",
       #    ref=REFS,
       #    wga=WGA,
       #    nowga=NOWGA
       #   #ext=["gz","gz.tbi"]
       #),
       ## pairwise gatk stats
       #expand(
       #    "results/{ref}/concordance/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
       #    ref=REFS,
       #    wga=WGA,
       #    nowga=NOWGA
       #   #ext=["gz","gz.tbi"]
       #),
        # final output from subset concordance test
       #"results/raw_subset_concord.gatk.csv"
        # begin main genotyping and imputation
       #expand(
       #   #"results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.{ext}",
       #   #"results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.{ext}",
       #   #"results/{ref}/vcf/{chrom}/powd_lowpass.imputed.snps.{chrom}.{ref}.vcf.{ext}",
       #   #"results/{ref}/vcf/imputed.list",
       #   #"results/{ref}/vcf/combine/powd_lowpass.imputed.snps.{ref}.vcf.{ext}",
       #    "results/{ref}/vcf/combine/powd_lowpass.dr208.imputed.snps.{ref}.vcf.{ext}",
       #    ref=REFS,
       #   #chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M CORRECT?
       #    ext=["gz","gz.tbi"]
       #),
       #expand(
       #    "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
       #    ref=REFS,
       #    wga=WGA,
       #    nowga=NOWGA
       #   #ext=["gz","gz.tbi"]
       #),
        "results/imputed_subset_concord.gatk.csv"
            

include: "rules/subset_concord.smk"
include: "rules/genotype.smk"
include: "rules/impute.smk"
include: "rules/impute_concord.smk"

