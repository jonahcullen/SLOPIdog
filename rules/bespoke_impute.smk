# GTDN WITH CHR25 USING NOT PHASED TARGET AND gl INSTEAD OF gt

rule split_phased_caravaggio_ref:
    output:
        ref_vcf = "results/{ref}/caravaggio_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        ref_tbi = "results/{ref}/caravaggio_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    params:
       #phased_pop = lambda wildcards, input: config['phased'][wildcards.ref]
       #phased_pop = "/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef/subset_ref.phased.vcf.gz"
        phased_pop = "/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef/caravaggio_ref.phased.vcf.gz"
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

rule beagle40_impute_gl:
    input:
       #phase_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz",
       #phase_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.phased.{chrom}.{ref}.vcf.gz.tbi",
        snps_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz",
        snps_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz.tbi",
       #ref_vcf   = "results/{ref}/caravaggio_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
       #ref_tbi   = "results/{ref}/caravaggio_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
        ref_vcf   = "results/{ref}/wakame_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz",
        ref_tbi   = "results/{ref}/wakame_pop/{chrom}/joint_genotype.{ref}.snps.{chrom}.vcf.gz.tbi"
    output:
        imputed_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.wakame_imputed.snps.{chrom}.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.wakame_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
    threads: 24
    resources:
        time   = 4320,
        mem_mb = 248000
    shell:
        '''
            java -Xms10g -Xmx240g -jar /opt/slopi/src/beagle4/beagle.r1399.jar \
                gl={input.snps_vcf} \
                ref={input.ref_vcf} \
                nthreads={threads} \
                impute=true \
                window=50000 \
                overlap=20000 \
                phase-its=12 \
                burnin-its=12 \
                out={params.prefix}

            tabix -p vcf {params.prefix}.vcf.gz
        '''

rule filter_imputed_wakame:
    input:
       #sorted_vcf = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz",
       #sorted_tbi = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz.tbi",
        imputed_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.wakame_imputed.snps.{chrom}.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.wakame_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.fltr_wakame_imputed.snps.{chrom}.{ref}.vcf.gz",
        filt_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.fltr_wakame_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools filter \
                -i 'INFO/DR2[*] >= 0.8' \
                -Oz -o {output.filt_vcf} \
                {input.imputed_vcf}

            tabix -p vcf {output.filt_vcf}
        '''

rule filter_imputed_caravaggio:
    input:
       #sorted_vcf = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz",
       #sorted_tbi = "results/{ref}/vcf/{breed}/combine/{breed}.cyrano_imputed.snps.{ref}.vcf.gz.tbi",
        imputed_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.caravaggio_imputed.snps.{chrom}.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.caravaggio_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/{breed}/{chrom}/{breed}.fltr_caravaggio_imputed.snps.{chrom}.{ref}.vcf.gz",
        filt_tbi = "results/{ref}/vcf/{breed}/{chrom}/{breed}.fltr_caravaggio_imputed.snps.{chrom}.{ref}.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools filter \
                -i 'INFO/DR2[*] >= 0.8' \
                -Oz -o {output.filt_vcf} \
                {input.imputed_vcf}

            tabix -p vcf {output.filt_vcf}
        '''

# GTDN with CHR25
rule split_phased_gtdn:
    output:
        ref_vcf = "results/{ref}/gtdn_pop/chr25/gtdn.{ref}.snps.chr25.vcf.gz",
        ref_tbi = "results/{ref}/gtdn_pop/chr25/gtdn.{ref}.snps.chr25.vcf.gz.tbi"
    params:
       #phased_pop = lambda wildcards, input: config['phased'][wildcards.ref]
        phased_pop = "/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef/gtdn/gtdn_ref.phased.vcf.gz"
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools view \
                -r chr25 \
                -Oz -o {output.ref_vcf} \
                {params.phased_pop}

            tabix -p vcf {output.ref_vcf}
        '''

rule beagle40_impute_gtdn_chr25:
    input:
        phase_vcf = "results/{ref}/vcf/gtdn/chr25/gtdn.snps.phased.chr25.{ref}.vcf.gz",
        phase_tbi = "results/{ref}/vcf/gtdn/chr25/gtdn.snps.phased.chr25.{ref}.vcf.gz.tbi",
        ref_vcf = "results/{ref}/gtdn_pop/chr25/gtdn.{ref}.snps.chr25.vcf.gz",
        ref_tbi = "results/{ref}/gtdn_pop/chr25/gtdn.{ref}.snps.chr25.vcf.gz.tbi"
    output:
        imputed_vcf = "results/{ref}/vcf/gtdn/chr25/gtdn.gtdn_ref_imputed.snps.chr25.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/gtdn/chr25/gtdn.gtdn_ref_imputed.snps.chr25.{ref}.vcf.gz.tbi",
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
    threads: 24
    resources:
        time   = 2160,
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

rule filter_imputed_gtdn:
    input:
        imputed_vcf = "results/{ref}/vcf/gtdn/chr25/gtdn.gtdn_ref_imputed.snps.chr25.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/gtdn/chr25/gtdn.gtdn_ref_imputed.snps.chr25.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/gtdn/filter/gtdn.fltr_gtdn_imputed_chr25.snps.{ref}.20231018.vcf.gz",
        filt_tbi = "results/{ref}/vcf/gtdn/filter/gtdn.fltr_gtdn_imputed_chr25.snps.{ref}.20231018.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools filter \
                -i 'INFO/DR2[*] >= 0.8' \
                -Oz -o {output.filt_vcf} \
                {input.imputed_vcf}

            tabix -p vcf {output.filt_vcf}
        '''

# GSHP with CHR24
rule split_phased_gshp:
    output:
        ref_vcf = "results/{ref}/gshp_pop/chr24/gshp.{ref}.snps.chr24.vcf.gz",
        ref_tbi = "results/{ref}/gshp_pop/chr24/gshp.{ref}.snps.chr24.vcf.gz.tbi"
    params:
       #phased_pop = lambda wildcards, input: config['phased'][wildcards.ref]
        phased_pop = "/panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef/gshp/gshp_ref.phased.vcf.gz"
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools view \
                -r chr24 \
                -Oz -o {output.ref_vcf} \
                {params.phased_pop}

            tabix -p vcf {output.ref_vcf}
        '''

rule beagle40_impute_gshp_chr24:
    input:
        phase_vcf = "results/{ref}/vcf/gshp/chr24/gshp.snps.phased.chr24.{ref}.vcf.gz",
        phase_tbi = "results/{ref}/vcf/gshp/chr24/gshp.snps.phased.chr24.{ref}.vcf.gz.tbi",
        ref_vcf = "results/{ref}/gshp_pop/chr24/gshp.{ref}.snps.chr24.vcf.gz",
        ref_tbi = "results/{ref}/gshp_pop/chr24/gshp.{ref}.snps.chr24.vcf.gz.tbi"
    output:
        imputed_vcf = "results/{ref}/vcf/gshp/chr24/gshp.gshp_ref_imputed.snps.chr24.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/gshp/chr24/gshp.gshp_ref_imputed.snps.chr24.{ref}.vcf.gz.tbi",
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit(".",2)[0]
    threads: 24
    resources:
        time   = 2160,
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

rule filter_imputed_gshp:
    input:
        imputed_vcf = "results/{ref}/vcf/gshp/chr24/gshp.gshp_ref_imputed.snps.chr24.{ref}.vcf.gz",
        imputed_tbi = "results/{ref}/vcf/gshp/chr24/gshp.gshp_ref_imputed.snps.chr24.{ref}.vcf.gz.tbi",
    output:
        filt_vcf = "results/{ref}/vcf/gshp/filter/gshp.fltr_gshp_imputed_chr24.snps.{ref}.20231018.vcf.gz",
        filt_tbi = "results/{ref}/vcf/gshp/filter/gshp.fltr_gshp_imputed_chr24.snps.{ref}.20231018.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 60,
        mem_mb = 24000
    shell:
        '''
            bcftools filter \
                -i 'INFO/DR2[*] >= 0.8' \
                -Oz -o {output.filt_vcf} \
                {input.imputed_vcf}

            tabix -p vcf {output.filt_vcf}
        '''
