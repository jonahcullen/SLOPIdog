import os
import itertools
import pandas as pd
from pathlib import Path

singularity: config['sif']

rule all:
    input:
       ## begin main genotyping and imputation
        expand(
           #"results/{ref}/impute/{breed}/combine/{breed}.dr2_fltr_imputed.snps.{ref}.{date}.vcf.gz.covtsf",
           #"{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/merseberg.{ref}.{chrom}.phased.vcf.gz.tbi",
           #"results/{ref}/target/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz.tbi",
            "results/{ref}/impute/{breed}/{chrom}/{breed}.fltr_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
            bucket=config['bucket'],
            ref="cf4",
            breed="gtdn",
           #breed=config['breeds'],
           #chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M OR Y
            chrom="chr25",
            date=config['date'],
        )

# NOTE the discrep between ref UU_Cfam for phasing and cf4 for genotyping/imputing
# NEED TO COMBINE THESE SENSIBILY BEFORE RE-PREPARING/PHASING THE WHOLE
# MERSEBERG PANEL
#include: "rules/phasing.merseberg.smk"
include: "rules/genotype.smk"
include: "rules/impute.smk"

