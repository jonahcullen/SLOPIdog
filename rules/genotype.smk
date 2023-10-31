
rule pile_dogs:
    input:
        bam_list = "/scratch.global/friedlab_LOW_PILEUP/SLOPIdog/breed_bams/scratch_bams.{breed}.list"
    output:
        lowpass_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf",
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
        lowpass_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf",
    output:
        lowpass_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz",
        lowpass_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz.tbi",
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
        lowpass_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz",
        lowpass_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.{chrom}.{ref}.vcf.gz.tbi",
    output:
        snps_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz",
        snps_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz.tbi",
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
                -Oz -o {output.snps_vcf} \
                {input.lowpass_vcf}

            tabix -p vcf {output.snps_vcf}
        '''

rule phase_x_chrom:
    input:
        snps_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz",
        snps_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz.tbi",
    output:
        phase_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz",
        phase_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz.tbi",
    params:
        link_map     = config["link_map"],
        eff_pop_size = 200,
        window       = 120,
        overlap      = 10,
        out_prefix   = lambda wildcards, output: output.phase_vcf.rsplit(".", 2)[0],
    threads: 24
    resources:
        time   = 140,
        mem_mb = 120000
    shell:
        '''
            java -jar -Xmx110g /opt/slopi/src/beagle5/beagle.22Jul22.46e.jar \
                gt={input.snps_vcf} \
                ne={params.eff_pop_size} \
                nthreads={threads} \
                map={params.link_map} \
                window={params.window} \
                overlap={params.overlap} \
                out={params.out_prefix}

            tabix -p vcf {output.phase_vcf}
        '''

