
rule pile_dogs:
    input:
        bam_list = config['bam_list']
    output:
        lowpass_vcf = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]['fasta']
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
        lowpass_vcf = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf",
    output:
        lowpass_vcf = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz",
        lowpass_tbi = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz.tbi",
    threads: 6
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        '''
            bgzip --threads {threads} -c {input.lowpass_vcf} > {output.lowpass_vcf}
            tabix -p vcf {output.lowpass_vcf}
        '''

rule filter_snps:
    input:
        lowpass_vcf = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz",
        lowpass_tbi = "results/{ref}/target/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz.tbi",
    output:
        snps_vcf = "results/{ref}/target/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz",
        snps_tbi = "results/{ref}/target/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz.tbi",
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
                -Oz -o {output.snps_vcf} \
                {input.lowpass_vcf}

            tabix -p vcf {output.snps_vcf}
        '''

# the pre-phasing strategy to split up the phasing/imputation will not work as
# in order to get "correct" results need to use gl field instead of gt in
# imputation - correct answer was verified using gtdns with chr25
