def parse_chunks(chunks_file):
    """Returns a dict of chunk_id:(input_region, output_region)"""
    chunks = {}
    with open(chunks_file) as f:
        for line in f:
            fields = line.split()
            chunks[fields[0]] = (fields[2], fields[3])
    return chunks

def get_chunk_ids(wildcards):
    """The chunks file is the output of chunk_reference named chunked_chrom, get it here and strip the first row as chunk names"""
    chunks_file = checkpoints.chunk_reference.get(**wildcards).output.chunked_chrom
    return list(parse_chunks(chunks_file).keys())

# need to convert our existing linkage map to one that fits the actual linkage map format expected by glimpse2
rule convert_link_map:
    output:
        chrom_link_map = "results/glimpse/{ref}/{chrom}/{ref}.{chrom}.link.map"
    params:
        linkage_map = lambda wildcards: config['link_map'][wildcards.ref]
    threads: 4
    resources:
        time = 60,
        mem_mb = 4000
    shell:
        '''
            python scripts/convert_link_map.py -i {params.linkage_map} -o {output.chrom_link_map} -r {wildcards.chrom}
        '''

checkpoint chunk_reference:
    input:
        chrom_link_map = "results/glimpse/{ref}/{chrom}/{ref}.{chrom}.link.map"
    output:
        chunked_chrom = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv"
    params:
        phased_pop = lambda wildcards: config['refgen'][wildcards.ref]['phased_pop']
    threads: 4
    resources:
        time = 120,
        mem_mb = 16000
    shell:
        '''
            GLIMPSE2_chunk \
                --input {params.phased_pop} \
                --region {wildcards.chrom} \
                --map {input.chrom_link_map} \
                --window-mb 4.0 \
                --buffer-mb 0.5 \
                --threads {threads} \
                --sequential \
                --output {output.chunked_chrom}.raw.tsv
            python scripts/fix_chunk_triple_overlap.py \
                -i {output.chunked_chrom}.raw.tsv \
                -o {output.chunked_chrom}
            rm {output.chunked_chrom}.raw.tsv
        '''

# we have to do a lot of funky magic here, getting the name that will actually be output and renaming it to something sensible, so that later rules can just use the chunk id
rule split_reference:
    input:
        chunks_file = "results/glimpse/{ref}/{chrom}/chunks.{ref}.{chrom}.tsv",
        chrom_link_map = "results/glimpse/{ref}/{chrom}/{ref}.{chrom}.link.map"
    output:
        binary_ref = "results/glimpse/{ref}/{chrom}/split/binary_ref.{ref}.{chrom}.{chunk_id}.bin"
    params:
        phased_pop = lambda wildcards: config['refgen'][wildcards.ref]['phased_pop'],
        input_region = lambda wildcards: parse_chunks(f"results/glimpse/{wildcards.ref}/{wildcards.chrom}/chunks.{wildcards.ref}.{wildcards.chrom}.tsv")[wildcards.chunk_id][0],
        output_region = lambda wildcards: parse_chunks(f"results/glimpse/{wildcards.ref}/{wildcards.chrom}/chunks.{wildcards.ref}.{wildcards.chrom}.tsv")[wildcards.chunk_id][1],
        prefix = "results/glimpse/{ref}/{chrom}/split/binary_ref.{ref}.{chrom}",
        glimpse_output = lambda wildcards: ("results/glimpse/{ref}/{chrom}/split/binary_ref.{ref}.{chrom}".format(**wildcards) + "_" + parse_chunks(f"results/glimpse/{wildcards.ref}/{wildcards.chrom}/chunks.{wildcards.ref}.{wildcards.chrom}.tsv")[wildcards.chunk_id][0].replace(":", "_").replace("-", "_") + ".bin")
    threads: 4
    resources:
        time = 120,
        mem_mb = 16000
    shell:
        '''
            GLIMPSE2_split_reference \
                --reference {params.phased_pop} \
                --map {input.chrom_link_map} \
                --input-region {params.input_region} \
                --output-region {params.output_region} \
                --threads {threads} \
                --output {params.prefix}
            mv {params.glimpse_output} {output.binary_ref}
        '''

# we don't need the genetic map because it's written into the binary reference file, same for input/output regions
rule impute_chunk:
    input:
        binary_ref = "results/glimpse/{ref}/{chrom}/split/binary_ref.{ref}.{chrom}.{chunk_id}.bin",
        bam_list = config['bam_list']
    output:
        imputed_chunk = temp("results/glimpse/{ref}/{chrom}/imputed/imputed.{ref}.{chrom}.{chunk_id}.bcf"),
        imputed_chunk_index = temp("results/glimpse/{ref}/{chrom}/imputed/imputed.{ref}.{chrom}.{chunk_id}.bcf.csi")
    params:
        input_region = lambda wildcards: parse_chunks(f"results/glimpse/{wildcards.ref}/{wildcards.chrom}/chunks.{wildcards.ref}.{wildcards.chrom}.tsv")[wildcards.chunk_id][0],
        output_region = lambda wildcards: parse_chunks(f"results/glimpse/{wildcards.ref}/{wildcards.chrom}/chunks.{wildcards.ref}.{wildcards.chrom}.tsv")[wildcards.chunk_id][1],
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]['fasta']
    threads: 16
    resources:
        time = 720,
        mem_mb = 128000
    shell:
        '''
            GLIMPSE2_phase \
                --bam-list {input.bam_list} \
                --reference {input.binary_ref} \
                --fasta {params.ref_fa} \
                --impute-reference-only-variants \
                --threads {threads} \
                --output {output.imputed_chunk}
        '''

rule ligate_single_chrom:
    input:
        imputed_chunks = lambda wildcards: expand("results/glimpse/{ref}/{chrom}/imputed/imputed.{ref}.{chrom}.{chunk_id}.bcf", ref=wildcards.ref, chrom=wildcards.chrom, chunk_id=get_chunk_ids(wildcards)),
        imputed_chunk_index = lambda wildcards: expand("results/glimpse/{ref}/{chrom}/imputed/imputed.{ref}.{chrom}.{chunk_id}.bcf.csi", ref=wildcards.ref, chrom=wildcards.chrom, chunk_id=get_chunk_ids(wildcards))
    output:
        ligated = "results/glimpse/{ref}/{chrom}/imputed.{ref}.{chrom}.vcf.gz",
        chunk_list = temp("results/glimpse/{ref}/{chrom}/imputed_chunk_list.{ref}.{chrom}.txt")
    threads: 4
    resources:
        time = 240,
        mem_mb = 32000
    shell:
        '''
            printf '%s\n' {input.imputed_chunks} > {output.chunk_list}
            GLIMPSE2_ligate \
                --input {output.chunk_list} \
                --output {output.ligated} \
                --threads {threads}
        '''

rule filter_glimpse_chroms:
    input:
        imputed_chrom = "results/glimpse/{ref}/{chrom}/imputed.{ref}.{chrom}.vcf.gz"
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

rule concat_unfiltered:
    input:
        imputed_chroms = lambda wildcards: expand("results/glimpse/{ref}/{chrom}/imputed.{ref}.{chrom}.vcf.gz", ref=wildcards.ref, chrom=config['refgen'][wildcards.ref]['chroms'])
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.unfiltered.glimpse_imputed.vcf.gz",
        index = "results/glimpse/{ref}/{ref}.unfiltered.glimpse_imputed.vcf.gz.tbi",
        file_list = temp("results/glimpse/{ref}/{ref}_unfiltered.txt")
    threads: 8
    resources:
        time = 1440,
        mem_mb = 8000
    shell:
        '''
            printf '%s\n' {input.imputed_chroms} > {output.file_list}
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {output.file_list}
            tabix {output.merged_vcf}
        '''

rule concat_filtered:
    input:
        filtered_chroms = lambda wildcards: expand("results/glimpse/{ref}/{chrom}/{chrom}.{ref}.imputed.filtered.vcf.gz", ref=wildcards.ref, chrom=config['refgen'][wildcards.ref]['chroms'])
    output:
        merged_vcf = "results/glimpse/{ref}/{ref}.filtered.glimpse_imputed.vcf.gz",
        index = "results/glimpse/{ref}/{ref}.filtered.glimpse_imputed.vcf.gz.tbi",
        file_list = temp("results/glimpse/{ref}/{ref}_filtered.txt")
    threads: 8
    resources:
        time = 1440,
        mem_mb = 8000
    shell:
        '''
            printf '%s\n' {input.filtered_chroms} > {output.file_list}
            bcftools concat \
                -Oz -o {output.merged_vcf} \
                -f {output.file_list}
            tabix {output.merged_vcf}
        '''
