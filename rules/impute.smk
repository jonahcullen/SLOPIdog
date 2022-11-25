
rule split_phased_cf3:
    output:
        phased_chrom = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        phased_tbi   = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    params:
        phased_pop = lambda wildcards, input: config['phased'][wildcards.ref]
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools view \
                -r {wildcards.chrom} \
                -Oz -o {output.phased_chrom} \
                {params.phased_pop}

            tabix -p vcf {output.phased_chrom}
        '''

rule beagle40_impute:
    input:
        snps_vcf_gz  = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz",
        snps_vcf_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.snps.{chrom}.{ref}.vcf.gz.tbi",
        phased_chrom = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        phased_tbi   = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    output:
        imputed_vcf = "results/{ref}/vcf/{chrom}/powd_lowpass.imputed.snps.{chrom}.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/{chrom}/powd_lowpass.imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
    threads: 24
    resources:
        time   = 2880,
        mem_mb = 248000
    shell:
        '''
            java -jar /opt/slopi/src/beagle.r1399.jar \
                gl={input.snps_vcf_gz} \
                ref={input.phased_chrom} \
                nthreads={threads} \
                impute=true \
                window=50000 \
                overlap=20000 \
                phase-its=12 \
                burnin-its=12 \
                out={params.prefix}

            tabix -p vcf {params.prefix}.vcf.gz
        '''

rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
            "results/{ref}/vcf/{chrom}/powd_lowpass.imputed.snps.{chrom}.{ref}.vcf.{ext}",
            ref=REFS,
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M CORRECT?
            ext=["gz","gz.tbi"]
        ))
    output:
        sorted_list = "results/{ref}/vcf/combine/imputed.list",
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        vcfs = filter(lambda vcf: vcf.endswith('.gz') and (wildcards.ref in vcf), input.imputed_vcfs)
        # sort chrom natural function
        # https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
        def natural_sort(l): 
            convert = lambda text: int(text) if text.isdigit() else text.lower()
            alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
            return sorted(l, key=alphanum_key)
        # sort input
        nat_sort = natural_sort(vcfs)
        # write to file
        with open(output.sorted_list,'w') as out:
            out.write('\n'.join(nat_sort))

rule combine_imputed:
    input:
        sorted_list = "results/{ref}/vcf/combine/imputed.list",
        imputed_vcfs = sorted(expand(
            "results/{ref}/vcf/{chrom}/powd_lowpass.imputed.snps.{chrom}.{ref}.vcf.{ext}",
            ref=REFS,
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M CORRECT?
            ext=["gz","gz.tbi"]
        ))
    output:
        sorted_vcf = "results/{ref}/vcf/combine/powd_lowpass.imputed.snps.{ref}.vcf.gz",
        sorted_tbi = "results/{ref}/vcf/combine/powd_lowpass.imputed.snps.{ref}.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.sorted_vcf} \
                -f {input.sorted_list}
            
            tabix -p vcf {output.sorted_vcf}
        '''

rule filter_imputed:
    input:
        sorted_vcf = "results/{ref}/vcf/combine/powd_lowpass.imputed.snps.{ref}.vcf.gz",
        sorted_tbi = "results/{ref}/vcf/combine/powd_lowpass.imputed.snps.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/combine/powd_lowpass.dr208.imputed.snps.{ref}.vcf.gz",
        filt_tbi = "results/{ref}/vcf/combine/powd_lowpass.dr208.imputed.snps.{ref}.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools filter \
                -i 'INFO/DR2[*] >= 0.8' \
                -Oz -o {output.filt_vcf} \
                {input.sorted_vcf}

            tabix -p vcf {output.filt_vcf}
        '''
    



