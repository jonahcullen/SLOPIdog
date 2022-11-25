
rule extract_imputed_vcfs:
    input:
        filt_vcf = "results/{ref}/vcf/combine/powd_lowpass.dr208.imputed.snps.{ref}.vcf.gz",
        filt_tbi = "results/{ref}/vcf/combine/powd_lowpass.dr208.imputed.snps.{ref}.vcf.gz.tbi",
    output:
        lowpass_sample     = "results/{ref}/vcf/impute/{sample}/{sample}_lowpass.{ref}.vcf.gz",
        lowpass_sample_tbi = "results/{ref}/vcf/impute/{sample}/{sample}_lowpass.{ref}.vcf.gz.tbi"
    threads: 4
    resources:
        time   = 720,
        mem_mb = 12000
    shell:
        '''
            bcftools view -Oz -c1 \
                -s {wildcards.sample} \
                {input.filt_vcf} \
                > {output.lowpass_sample}

            tabix -p vcf {output.lowpass_sample}
        '''

rule pairwise_imputed_gatk:
    input:
        unpack(get_pairs)
    output:
        contin_met  = "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
        detail_met  = "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_detail_metrics",
        summary_met = "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_summary_metrics",
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

rule gather_imputed_gatk_summaries:
    input:
        summaries = expand(
            "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_summary_metrics",
            ref=["cf3","cf4"],
            wga=WGA,
            nowga=NOWGA
        ),
        contins = expand(
            "results/{ref}/concordance/imputed/gatk/{wga}_x_{nowga}_lowpass.genotype_concordance_contingency_metrics",
            ref=["cf3","cf4"],
            wga=WGA,
            nowga=NOWGA
        ),
    output:
        all_concords = "results/imputed_subset_concord.gatk.csv"
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
                    
