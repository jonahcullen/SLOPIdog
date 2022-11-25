
rule pairwise:
    input:
        lowpass_sample = "results/{ref}/vcf/{sample}/{sample}_lowpass.{ref}.vcf.gz"
    output:
        lowpass_comp = "results/{ref}/concordance/{wga}_{nowga}_lowpass.{ref}.stats.txt",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            echo HELLO
        '''
