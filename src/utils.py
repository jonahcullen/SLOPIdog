
def get_pairs(wildcards):
    return {
        'wga': f"results/{wildcards.ref}/vcf/impute/{wildcards.wga}/{wildcards.wga}_lowpass.{wildcards.ref}.vcf.gz", 
        'nowga': f"results/{wildcards.ref}/vcf/impute/{wildcards.nowga}/{wildcards.nowga}_lowpass.{wildcards.ref}.vcf.gz"
    }
