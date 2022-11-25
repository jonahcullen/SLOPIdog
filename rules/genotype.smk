
rule pile_dogs:
    input:
        bam_list = "/scratch.global/friedlab_BCFJOINT/lowpass_bams.{ref}.list"
    output:
        lowpass_vcf = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools mpileup \
                -Ou \
                -f {params.ref_fa} \
                -b {input.bam_list} \
                -r {wildcards.chrom} \
            | bcftools call -Oz -mv \
            | bcftools filter -s LowQual -e '%QUAL<20' > {output.lowpass_vcf}
        '''

rule zip_and_index:
    input:
        lowpass_vcf = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf",
    output:
        lowpass_vcf_gz  = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf.gz",
        lowpass_vcf_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf.gz.tbi",
    threads: 6
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        '''
            bgzip --threads {threads} -c {input.lowpass_vcf} > {output.lowpass_vcf_gz}
            tabix -p vcf {output.lowpass_vcf_gz}
        '''

rule filter_snps:
    input:
        lowpass_vcf_gz  = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf.gz",
        lowpass_vcf_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.{chrom}.{ref}.vcf.gz.tbi",
    output:
        snps_vcf_gz  = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz",
        snps_vcf_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz.tbi",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools view \
                -m 2 \
                -M 2 \
                -v snps \
                -r {wildcards.chrom} \
                -e ' GT="." ' \
                -Oz -o {output.snps_vcf_gz} \
                {input.lowpass_vcf_gz}

            tabix -p vcf {output.snps_vcf_gz}
        '''

#rule zip_and_index:
#    input:
#        snps_vcf = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf",
#    output:
#        snps_vcf_gz  = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz",
#        snps_vcf_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz.tbi",
#    threads: 6
#    resources:
#        time   = 720,
#        mem_mb = 12000
#    shell:
#        '''
#            bgzip --threads {threads} -c {input.snps_vcf} > {output.snps_vcf_gz}
#            tabix -p vcf {output.snps_vcf_gz}
#        '''
