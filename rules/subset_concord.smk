
rule pile_dogs_subset:
    input:
        bam_list = "/scratch.global/friedlab_BCFJOINT/concord_bams.{ref}.list"
    output:
        lowpass_vcf = "results/{ref}/vcf/subset_lowpass.{ref}.vcf",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]
    threads: 12
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            bcftools mpileup \
                -Ou \
                -f {params.ref_fa} \
                -b {input.bam_list} \
            | bcftools call -Oz -mv \
            | bcftools filter -s LowQual -e '%QUAL<20' > {output.lowpass_vcf}
        '''

rule zip_and_index_subset:
    input:
        lowpass_vcf = "results/{ref}/vcf/subset_lowpass.{ref}.vcf",
    output:
        lowpass_vcf_gz  = "results/{ref}/vcf/subset_lowpass.{ref}.vcf.gz",
        lowpass_vcf_tbi = "results/{ref}/vcf/subset_lowpass.{ref}.vcf.gz.tbi",
    threads: 12
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        '''
            bgzip -c {input.lowpass_vcf} > {output.lowpass_vcf_gz}
            tabix -p vcf {output.lowpass_vcf_gz}
        '''

rule extract_sample_vcfs:
    input:
        lowpass_vcf_gz  = "results/{ref}/vcf/subset_lowpass.{ref}.vcf.gz",
        lowpass_vcf_tbi = "results/{ref}/vcf/subset_lowpass.{ref}.vcf.gz.tbi",
    output:
        lowpass_sample     = "results/{ref}/vcf/{sample}/{sample}_lowpass.{ref}.vcf.gz",
        lowpass_sample_tbi = "results/{ref}/vcf/{sample}/{sample}_lowpass.{ref}.vcf.gz.tbi"
    threads: 4
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        '''
            bcftools view -Oz -c1 \
                -s {wildcards.sample} \
                {input.lowpass_vcf_gz} \
                > {output.lowpass_sample}

            tabix -p vcf {output.lowpass_sample}
        '''

rule pairwise_bcftools:
    input:
        unpack(get_pairs)
    output:
        lowpass_comp = "results/{ref}/concordance/raw/bcftools/{wga}_x_{nowga}_lowpass.{ref}.stats.txt",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref]
    threads: 4
    resources:
        time   = 240,
        mem_mb = 12000
    shell:
        '''
            bcftools stats {input.wga} {input.nowga} > {output.lowpass_comp}
        '''

rule pairwise_gatk:
    input:
        unpack(get_pairs)
    output:
        contin_met  = "results/{ref}/concordance/raw/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
        detail_met  = "results/{ref}/concordance/raw/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_detail_metrics",
        summary_met = "results/{ref}/concordance/raw/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_summary_metrics",
    params:
        ref_fa = lambda wildcards, input: config['refgen'][wildcards.ref],
        base_out = lambda wildcards, output: output.contin_met.rsplit(".",1)[0]
    threads: 4
    resources:
        time   = 240,
        mem_mb = 12000
    shell:
        '''
            gatk GenotypeConcordance \
                --CALL_VCF {input.wga} \
                --CALL_SAMPLE {wildcards.wga} \
                --OUTPUT {params.base_out} \
                --TRUTH_VCF {input.nowga} \
                --TRUTH_SAMPLE {wildcards.nowga}
        '''

rule gather_gatk_summaries:
    input:
        summaries = expand(
            "results/{ref}/concordance/raw/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_summary_metrics",
            ref=["cf3","cf4"],
            wga=WGA,
            nowga=NOWGA
        ),
        contins = expand(
            "results/{ref}/concordance/raw/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
            ref=["cf3","cf4"],
            wga=WGA,
            nowga=NOWGA
        ),
    output:
        all_concords = "results/raw_subset_concord.gatk.csv"
    run:
        from collections import defaultdict
        # read in summary metrics and parse snp/indel concordances
        d = defaultdict(list)
        for s in input.summaries:
            with open(s,'r') as f:
                for line in f:
                    # get ref from file name
                    ref = s.split('/')[1]
                    # grab snp concordance
                    if line.startswith("SNP"):
                        line = line.strip().split('\t')
                        # check if samples same
                        same = "no"
                        if line[1].split("_")[0] == line[2].split("_")[0]:
                            same = "yes"
                        key = f"SNP_{ref}_{line[1]}_{line[2]}"
                        d[key].extend([line[1],line[2],same,ref,"SNP",line[-2],line[-1]])
                    # grab indel concordance
                    elif line.startswith("INDEL"):
                        line = line.strip().split('\t')
                        key = f"INDEL_{ref}_{line[1]}_{line[2]}"
                        d[key].extend([line[1],line[2],same,ref,"INDEL",line[-2],line[-1]])

        # read in contingency metrics and parse snp/indel counts
        for c in input.contins:
            with open(c,'r') as f:
                for line in f:
                    # get ref from file name
                    ref = c.split('/')[1]
                    # grab snp counts
                    if line.startswith("SNP"):
                        line = line.strip().split('\t')
                        key = f"SNP_{ref}_{line[1]}_{line[2]}"
                        var_ct = sum(map(int,line[3:]))
                        d[key].append(var_ct)
                    elif line.startswith("INDEL"):
                        line = line.strip().split('\t')
                        key = f"INDEL_{ref}_{line[1]}_{line[2]}"
                        var_ct = sum(map(int,line[3:]))
                        d[key].append(var_ct)
        
        # write d to output
        with open(output.all_concords,'w') as out:
            print("truth_sample,call_sample,same,ref,type,geno_concord,non_ref_geno_concord,count",file=out)
            for v in d.values():
                print(','.join(str(i) for i in v),file=out)
                    
