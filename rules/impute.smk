
rule split_phased_ref:
    output:
        ref_vcf = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        ref_tbi = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    params:
       #phased_pop = lambda wildcards, input: config['phased'][wildcards.ref]
        phased_pop = "/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef/subset_ref.phased.vcf.gz"
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools view \
                -r {wildcards.chrom} \
                -Oz -o {output.ref_vcf} \
                {params.phased_pop}

            tabix -p vcf {output.ref_vcf}
        '''

#rule beagle54_impute:
#    input:
#       #snps_vcf  = "results/{ref}/vcf/{chrom}/all_lowpass.snps.{chrom}.{ref}.vcf.gz",
#       #snps_vcf = "results/{ref}/vcf/{chrom}/all_lowpass.snps.{chrom}.{ref}.vcf.gz.tbi",
#        phase_vcf = "results/{ref}/vcf/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz",
#        phase_tbi = "results/{ref}/vcf/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz.tbi",
#        ref_vcf   = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
#        ref_tbi   = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
#       #phase_vcf    = "results/{ref}/vcf/{chrom}/all_lowpass.snps.phased.{chrom}.{ref}.vcf.gz",
#       #phase_tbi    = "results/{ref}/vcf/{chrom}/all_lowpass.snps.phased.{chrom}.{ref}.vcf.gz.tbi",
#       #phased_chrom = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
#       #phased_tbi   = "results/{ref}/pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
#    output:
#       #imputed_vcf = "results/{ref}/vcf/{chrom}/all_lowpass.TEST_imputed.snps.{chrom}.{ref}.vcf.gz",
#       #imputed_tbi = "results/{ref}/vcf/{chrom}/all_lowpass.TEST_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
#        imputed_vcf = "results/{ref}/vcf/{chrom}/{breed}.cyrano_beagle5_burn12_imputed.snps.{chrom}.{ref}.vcf.gz",
#        imputed_tbi = "results/{ref}/vcf/{chrom}/{breed}.cyrano_beagle5_burn12_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
#    params:
#        link_map = config["link_map"],
#        prefix   = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
#    threads: 24
#    resources:
#        time   = 2880,
#        mem_mb = 248000
#    shell:
#        '''
#            java -Xms10g -Xmx220g -jar /opt/slopi/src/beagle5/beagle.22Jul22.46e.jar \
#                gt={input.phase_vcf} \
#                ref={input.phased_chrom} \
#                map={params.link_map} \
#                nthreads={threads} \
#                impute=true \
#                window=40.0 \
#                overlap=2.0 \
#                burnin=12 \
#                iterations=12 \
#                out={params.prefix}
#
#            tabix -p vcf {params.prefix}.vcf.gz
#        '''

rule beagle40_impute:
    input:
        phase_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz",
        phase_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz.tbi",
        ref_vcf   = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        ref_tbi   = "results/{ref}/cyrano_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    output:
        imputed_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.cyrano_imputed.snps.{chrom}.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.cyrano_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
    threads: 24
    resources:
        time   = 4320,
        mem_mb = 248000
    shell:
        '''
            java -Xms10g -Xmx240g -jar /opt/slopi/src/beagle4/beagle.r1399.jar \
                gt={input.phase_vcf} \
                ref={input.ref_vcf} \
                nthreads={threads} \
                impute=true \
                window=50000 \
                overlap=20000 \
                usephase=true \
                phase-its=0 \
                burnin-its=0 \
                out={params.prefix}

            tabix -p vcf {params.prefix}.vcf.gz
        '''

rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
            "results/{ref}/vcf/{breed}/{chrom}/{breed}.cyrano_imputed.snps.{chrom}.{ref}.vcf.{ext}",
            ref=REFS,
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"], # NO CHROM M OR Y
            ext=["gz","gz.tbi"],
            breed=config["breeds"]
        ))
    output:
       #sorted_list = "results/{ref}/vcf/combine/imputed.list",
        sorted_list = "results/{ref}/vcf/{breed}/combine/{breed}_imputed.list",
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        vcfs = filter(lambda vcf: vcf.endswith('.gz') and (wildcards.ref in vcf) and (wildcards.breed in vcf), input.imputed_vcfs)
        print(vcfs)
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
        sorted_list = "results/{ref}/vcf/{breed}/combine/{breed}_imputed.list",
    output:
        sorted_vcf = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz",
        sorted_tbi = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 480,
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
        sorted_vcf = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz",
        sorted_tbi = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/{breed}/combine/{breed}.fltr_cyrano_imputed.snps.{ref}.20231018.vcf.gz",
        filt_tbi = "results/{ref}/vcf/{breed}/combine/{breed}.fltr_cyrano_imputed.snps.{ref}.20231018.vcf.gz.tbi",
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
    
