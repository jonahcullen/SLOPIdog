# THIS IS TECHNICALLY NOT ALLOWED
# NEED TO REWRITE THIS TO NOT REQUIRE REDISTRIBUTING STUFF THAT ISN'T OURS, BUT IS OKAY* FOR AN INTERNAL TOOL
# *not ok

# takes no inputs, as chunks are determined from the phased reference (should that be an input instead of a param???)
rule chunk_reference:
    output:
        chunked_chrom = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv"
    params:
        phased_pop = lambda wildcards, input: config['refgen'][wildcards.ref]['phased']
    threads: 4
    resources:
        time: 120,
        mem_mb: 16000
    shell:
        '''
            GLIMPSE_static/GLIMPSE_chunk_static --input {params.phased_pop} --region {wildcards.chrom} --window-size 1000000 --window-count 1000 --buffer-size 250000 --buffer-count 250 --output {output.chunked_chrom}
        '''


rule impute_chunk:
    input:
        # previously generated chromosome chunk file
        chunked_chrom = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv"
        # genotype.smk output panel to be imputed
        target_panel = "results/{ref}/target/{breed}/{chrom}/{breed}.snps.{chrom}.{ref}.vcf.gz"
    output:
        chrom_manifest = "results/glimpse/{ref}/{chrom}/{chrom}_manifest.txt"
    params:
        # healthy phased population
        phased_pop = lambda wildcards, input: config['refgen'][wildcards.ref]['phased']
        # config-controlled genetic map
        linkage_map = lambda wildcards, input: config['link_map'][wildcards.ref]
    threads: 16
    resources:
        time: 2880
        mem_mb: 48000
    # iterate over each part of the chunk file and run imputation, outputting a file for each, then output a manifest of generated files
    shell:
        '''
            ITER=1
            while IFS="" read -r LINE || [ -n "$LINE" ];
            do
                printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
                IRG=$(echo $LINE | cut -d" " -f3)
                ORG=$(echo $LINE | cut -d" " -f4)
                OUT=results/glimpse/{wildcards.ref}/{wildcards.chrom}/imputed_chunks/{wildcards.chrom}.imputed.${ITER}.bcf
                GLIMPSE_static/GLIMPSE_phase_static --input {input.target_panel} --reference {params.phased_pop} --map {params.linkage_map} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
                bcftools index -f ${OUT}
                ITER=$(($ITER + 1))
            done < {input.chunked_chrom}
            ls results/glimpse/{wildcards.ref}/{wildcards.chrom}/imputed_chunks/*.bcf > {output.chrom_manifest}
        '''

rule ligate_chroms:
    input:
        chrom_manifest = "results/glimpse/{ref}/{chrom}/{chrom}_manifest.txt"
    output:
        imputed_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf"
        index = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf.csi"
    threads: 4
    resources:
        time: 480
        mem_mb: 64000
    shell:
        '''
            GLIMPSE_static/GLIMPSE_ligate_static --input {input.chrom_manifest} --output {output.imputed_chrom}
            bcftools index -f {output.imputed_chrom}
        '''

# THIS RULE DOESN'T ACTUALLY SORT ANYTHING CURRENTLY
rule sort_chroms:
    input:
        # expand to the list of all chroms we have
        imputed_chroms = expand("results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf", 
            ref = config['refs'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"],
            breed=config["breeds"])
    output:
        sorted_chrom_bcfs = "results/glimpse/{ref}/{ref}_sorted_bcfs.txt"
    threads: 4
    resources:
        time: 120
        mem_mb: 4000
    run:
        # yoinked from impute.smk
        # make sure it ends with the correct extension and that we're only dealing with same-ref files
        bcfs = filter(lambda bcf: bcf.endswith('.bcf') and (wildcards.ref in bcf), input.imputed_chroms)
        with open({output.sorted_chrom_bcfs}, "wt") as outfile:
            for bcf in bcfs:
                # get the region name
                chrom = bcf.split('/')
                chrom = chrom[len(pro_bcf) - 1].split('.')
                chrom = chrom[0]
                print(bcf, file = outfile)

rule concat_chroms:
    input:
        bcf_list = "results/glimpse/{ref}/{ref}_sorted_bcfs.txt"
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.glimpse_imputed.vcf.gz"
        index = "results/glimpse/{ref}/{ref}.glimpse_imputed.vcf.gz.tbi"
    threads: 8
    resources:
        time: 1440
        mem_mb: 64000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {input.bcf_list}
            tabix {output.merged_vcf}
        '''
                
rule filter_glimpse_imputed:
    input:
        merged_vcf = "results/glimpse/{ref}/{ref}.glimpse_imputed.vcf.gz"
    output:
        filtered_vcf = "results/glimpse/{ref}/{ref}.glimpse_imputed.filtered.vcf.gz"
        index = "results/glimpse/{ref}/{ref}.glimpse_imputed.filtered.vcf.gz.tbi"
    threads: 4
    resources:
        time: 720
        mem_mb: 24000
    shell:
        '''
            bcftools filter -e 'INFO/INFO<0.8' -Oz -o {output.filtered_vcf} {input.merged_vcf}
            tabix {output.filtered_vcf}
        '''

rule phase_filtered_glimpse:
    input:
        filtered_vcf = "results/glimpse/{ref}/{ref}.glimpse_imputed.filtered.vcf.gz"
    output:
        phased_vcf = "results/glimpse/{ref}/{ref}.glimpse_imputed.filtered.phased.vcf.gz"
        index = "results/glimpse/{ref}/{ref}.glimpse_imputed.filtered.phased.vcf.gz.tbi"
        log = "results/glimpse/{ref}/{ref}.phasing_log.txt"
    threads: 8
    resources:
        time: 1440
        mem_mb: 64000
    shell:
        '''
            GLIMPSE_static/GLIMPSE_sample_static --input {input.filtered_vcf} --solve --output {output.phased_vcf} --log {output.log}
            tabix {output.phased_vcf}
        '''