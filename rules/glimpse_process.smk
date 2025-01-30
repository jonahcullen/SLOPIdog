# THIS IS TECHNICALLY NOT ALLOWED
# NEED TO REWRITE THIS TO NOT REQUIRE REDISTRIBUTING STUFF THAT ISN'T OURS, BUT IS OKAY* FOR AN INTERNAL TOOL
# *not ok

# takes no inputs, as chunks are determined from the phased reference (should that be an input instead of a param???)
rule chunk_reference:
    output:
        chunked_chrom = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv"
    params:
        phased_pop = lambda wildcards, input: config['refgen'][wildcards.ref]['phased_pop']
    threads: 4
    resources:
        time = 120,
        mem_mb = 16000
    shell:
        '''
            mkdir results/glimpse/{wildcards.ref}/{wildcards.chrom}/imputed_chunks
            GLIMPSE_static/GLIMPSE_chunk_static --input {params.phased_pop} --region {wildcards.chrom} --window-size 1000000 --window-count 1000 --buffer-size 250000 --buffer-count 250 --output {output.chunked_chrom}
        '''

rule convert_link_map:
    output:
        chrom_link_map = "results/glimpse/{ref}/{chrom}/{ref}.{chrom}.link.map"
    params:
        linkage_map = lambda wildcards, input: config['link_map'][wildcards.ref]
    threads: 4
    resources:
        time = 60,
        mem_mb = 4000
    shell:
        '''
            python GLIMPSE_static/convert_link_map.py -i {params.linkage_map} -o {output.chrom_link_map} -r {wildcards.chrom}
        '''

rule impute_chunk:
    input:
        # previously generated chromosome chunk file
        chunked_chrom = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv",
        # genotype.smk output panel to be imputed
        target_panel = "results/{ref}/target/all/{chrom}/all.snps.{chrom}.{ref}.vcf.gz",
        linkage_map = "results/glimpse/{ref}/{chrom}/{ref}.{chrom}.link.map"
    output:
        chrom_manifest = "results/glimpse/{ref}/{chrom}/{chrom}_manifest.txt"
    params:
        # healthy phased population
        phased_pop = lambda wildcards, input: config['refgen'][wildcards.ref]['phased_pop']
    threads: 32
    resources:
        time = 720,
        mem_mb = 48000
    # iterate over each part of the chunk file and run imputation, outputting a file for each, then output a manifest of generated files
    shell:
        '''
            ITER=1
            while IFS="" read -r LINE || [ -n "$LINE" ];
            do
                printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
                IRG=$(echo $LINE | cut -d" " -f3)
                ORG=$(echo $LINE | cut -d" " -f4)
                OUT=results/glimpse/{wildcards.ref}/{wildcards.chrom}/imputed_chunks/{wildcards.chrom}.imputed.$ITER.bcf
                GLIMPSE_static/GLIMPSE_phase_static --input {input.target_panel} --reference {params.phased_pop} --map {input.linkage_map} --input-region $IRG --output-region $ORG --output $OUT --thread 32
                bcftools index -f $OUT
                ITER=$(($ITER + 1))
            done < {input.chunked_chrom}
            ls results/glimpse/{wildcards.ref}/{wildcards.chrom}/imputed_chunks/*.bcf > {output.chrom_manifest}
        '''

rule ligate_chrom:
    input:
        chrom_manifest = "results/glimpse/{ref}/{chrom}/{chrom}_manifest.txt"
    output:
        imputed_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf",
        index = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf.csi"
    threads: 4
    resources:
        time = 480,
        mem_mb = 64000
    shell:
        '''
            GLIMPSE_static/GLIMPSE_ligate_static --input {input.chrom_manifest} --output {output.imputed_chrom}
            bcftools index -f {output.imputed_chrom}
        '''

rule filter_glimpse_chroms:
    input:
        imputed_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf"
    output:
        filtered_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.filtered.vcf.gz",
        index = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.filtered.vcf.gz.tbi"
    threads: 4
    resources:
        time = 720,
        mem_mb = 24000
    shell:
        '''
            bcftools filter -e 'INFO/INFO<0.8' -Oz -o {output.filtered_chrom} {input.imputed_chrom}
            tabix {output.filtered_chrom}
        '''

rule phase_filtered_glimpse_chroms:
    input:
        filtered_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.filtered.vcf.gz"
    output:
        phased_chrom = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.phased.filtered.vcf.gz",
        index = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.phased.filtered.vcf.gz.tbi",
        log = "results/glimpse/{ref}/{chrom}/{chrom}.{ref}.phasing_log.txt"
    threads: 8
    resources:
        time = 1440,
        mem_mb = 64000
    shell:
        '''
            GLIMPSE_static/GLIMPSE_sample_static --input {input.filtered_chrom} --solve --output {output.phased_chrom} --log {output.log}
            tabix {output.phased_chrom}
        '''

# this is certainly adaptable to be less obtuse (making a separate manifest for all three deliverables is a... choice)
rule sort_chroms:
    input:
        # expand to the list of all chroms we have
        imputed_chroms = expand("results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.bcf", 
            ref = config['refs'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"]),
        filtered_chroms = expand("results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.filtered.vcf.gz", 
            ref = config['refs'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"]),
        phased_chroms = expand("results/glimpse/{ref}/{chrom}/{chrom}.{ref}.phased.filtered.vcf.gz", 
            ref = config['refs'],
            chrom=[f"chr{i}" for i in range(1,39)] + ["chrX"])
    output:
        sorted_unfiltered = "results/glimpse/{ref}/{ref}_unfiltered.txt",
        sorted_filtered = "results/glimpse/{ref}/{ref}_filtered.txt",
        sorted_phased = "results/glimpse/{ref}/{ref}_phased.txt"
    threads: 4
    resources:
        time = 120,
        mem_mb = 4000
    run:
        # yoinked from impute.smk, thanks Jonah!
        # drop indices from input (do we need to do this?)
        unfiltered_bcfs = filter(lambda vcf: vcf.endswith('.bcf') and (wildcards.ref in vcf), input.imputed_chroms)
        print(unfiltered_bcfs)
        filtered_vcfs = filter(lambda vcf: vcf.endswith('.gz') and (wildcards.ref in vcf), input.filtered_chroms)
        print(filtered_vcfs)
        phased_vcfs = filter(lambda vcf: vcf.endswith('.gz') and (wildcards.ref in vcf), input.phased_chroms)
        print(phased_vcfs)
        # sort chrom natural function (ooooh cool!)
        # https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
        def natural_sort(l): 
            convert = lambda text: int(text) if text.isdigit() else text.lower()
            alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
            return sorted(l, key=alphanum_key)
        # sort input
        nat_sort_unfiltered = natural_sort(unfiltered_bcfs)
        nat_sort_filtered = natural_sort(filtered_vcfs)
        nat_sort_phased = natural_sort(phased_vcfs)
        # write to file
        with open(output.sorted_unfiltered,'w') as out:
            out.write('\n'.join(nat_sort_unfiltered))
        with open(output.sorted_filtered,'w') as out:
            out.write('\n'.join(nat_sort_filtered))
        with open(output.sorted_phased,'w') as out:
            out.write('\n'.join(nat_sort_phased))


# very inelegant because my elegant one failed, these should be one rule probably generated above with "stage"
rule concat_unfiltered:
    input:
        file_list = "results/glimpse/{ref}/{ref}_unfiltered.txt"
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.unfiltered.glimpse_imputed.vcf.gz",
        index = "results/glimpse/{ref}/{ref}.unfiltered.glimpse_imputed.vcf.gz.tbi"
    threads: 8
    resources:
        time = 1440,
        mem_mb = 64000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {input.file_list}
            tabix {output.merged_vcf}
        '''

rule concat_filtered:
    input:
        file_list = "results/glimpse/{ref}/{ref}_filtered.txt"
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.filtered.glimpse_imputed.vcf.gz",
        index = "results/glimpse/{ref}/{ref}.filtered.glimpse_imputed.vcf.gz.tbi"
    threads: 8
    resources:
        time = 1440,
        mem_mb = 64000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {input.file_list}
            tabix {output.merged_vcf}
        '''

rule concat_phased:
    input:
        file_list = "results/glimpse/{ref}/{ref}_phased.txt"
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.phased.glimpse_imputed.vcf.gz",
        index = "results/glimpse/{ref}/{ref}.phased.glimpse_imputed.vcf.gz.tbi"
    threads: 8
    resources:
        time = 1440,
        mem_mb = 64000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {input.file_list}
            tabix {output.merged_vcf}
        '''