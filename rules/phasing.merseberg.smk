
# subset the joint call vcf from 202308 to charlemagne using specific dog
# breeds plus 17 others - see /panfs/jay/groups/0/fried255/shared/gatk4_workflow/LowPass/SlimRef
# for jobs and sample lists. then used bcftools to recalculate AC/AF
# bcftools +fill-tags /scratch.global/friedlab_EALDBERHT/charlemagne_ref.vcf.gz -Oz -o /scratch.global/friedlab_EALDBERHT/FILL_TAGS_charlemagne_ref.vcf.gz
# Simultaneously created a workflow to use vcffixup from vcflib
# see /panfs/jay/groups/0/fried255/cull0084/projects/Misc/SnakeLines/RecalcAF

# NEED to add a rule for bcftools +fill-tags as the above steps were completed
# prior to starting pipeline
rule fix_ploidy_x_chrom:
    input:
       #final_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/final_gather/merseberg.af_fix.vcf.gz",
        final_vcf = "/scratch.global/friedlab_LOWPASS/merseberg.UU_Cfam_GSD_1.0_ROSY.chr25.fill_tags.vcf.gz",
    output:
        tmp_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/force_fix_ploidy/merseberg.{ref}.{chrom}.vcf.gz",
        tmp_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/force_fix_ploidy/merseberg.{ref}.{chrom}.vcf.gz.tbi",
    threads: 4
    resources:
         time   = 2880,
         mem_mb = 16000
    shell:
        '''
            # change chr39 to X to extract
            chrom={wildcards.chrom}
            if [ $chrom = "chr39" ]; then
                chrom=chrX
            fi
            
            bcftools +fixploidy \
                --regions $chrom \
                -Oz \
                -o {output.tmp_vcf} \
                {input.final_vcf} \
                -- \
                --force-ploidy 2

            tabix -p vcf {output.tmp_vcf}
        '''

rule filter_x_chrom:
    input:
        tmp_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/force_fix_ploidy/merseberg.{ref}.{chrom}.vcf.gz",
        tmp_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/force_fix_ploidy/merseberg.{ref}.{chrom}.vcf.gz.tbi",
    output:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/merseberg.{ref}.{chrom}.vcf.gz",
    threads: 4
    resources:
         time   = 60,
         mem_mb = 16000
    shell:
        '''
            bcftools view \
                --min-alleles 2 \
                --max-alleles 2 \
                --types snps \
                --include 'INFO/AF > 0.001071811 && GT!="."' \
                {input.tmp_vcf} \
            | \
            bcftools filter \
                --exclude "F_MISSING > 0.01 || FILTER!='PASS'" \
                -Oz \
                -o {output.chrom_vcf}
        '''
                #--exclude "F_MISSING > 0.01 || FILTER!='PASS' || FILTER!='ExcessHet' || FILTER!='VQSRTrancheSNP99.00to99.50' || FILTER!='VQSRTrancheSNP99.50to99.90'" \
           #bcftools view \
           #    --min-alleles 2 \
           #    --max-alleles 2 \
           #    --types snps \
           #    --exclude ' GT="." ' \
           #    --include 'INFO/AF > 0.00134770889487871' \
           #    {input.tmp_vcf} \
           #| \
           #bcftools filter \
           #    --exclude "F_MISSING > 0.01" \
           #    -f PASS \
           #    -Oz \
           #    -o {output.chrom_vcf}

rule phase_x_chrom:
    input:
        chrom_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/filtered/merseberg.{ref}.{chrom}.vcf.gz",
    output:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/merseberg.{ref}.{chrom}.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/merseberg.{ref}.{chrom}.phased.vcf.gz.tbi",
    params:
        link_map     = "/home/refgen/dog/canfam4/canFam4.linkage.map",
        eff_pop_size = 200,
        window       = 120,
        overlap      = 10,
        out_prefix   = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/merseberg.{ref}.{chrom}.phased",
    threads: 24
    resources:
        time   = 1440,
        mem_mb = 248000
    shell:
        '''
            java -jar -Xmx246g /opt/slopi/src/beagle5/beagle.22Jul22.46e.jar \
                gt={input.chrom_vcf} \
                ne={params.eff_pop_size} \
                nthreads={threads} \
                map={params.link_map} \
                window={params.window} \
                overlap={params.overlap} \
                out={params.out_prefix}

            tabix -p vcf {output.phase_vcf}
        '''

rule concat_phased:
    input:
        phase_vcf = S3.remote(expand(
            "{bucket}/wgs/pipeline/{ref}/{date}/phasing/phased/charlemagne.{ref}.{chrom}.phased.vcf.gz",
            bucket=config['bucket'],
            ref=config['ref'],
            date=config['date'],
            chrom=[f"chr{i}" for i in range(1,39+1)]
        )),
    output:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/charlemagne.{ref}.{date}.snps.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/charlemagne.{ref}.{date}.snps.phased.vcf.gz.tbi",
    threads: 4
    resources:
        time   = 720,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
                -Oz -o {output.phase_vcf} \
                {input.phase_vcf}

            tabix -p vcf {output.phase_vcf}
        '''

rule gbindex_phase_vcf:
    input:
        phase_vcf = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/charlemagne.{ref}.{date}.snps.phased.vcf.gz",
        phase_tbi = "{bucket}/wgs/pipeline/{ref}/{date}/phasing/charlemagne.{ref}.{date}.snps.phased.vcf.gz.tbi",
    output:
        "{bucket}/wgs/pipeline/{ref}/{date}/phasing/charlemagne.{ref}.{date}.snps.phased.vcf.gz.covtsf",
    params:
        ref_dir = os.path.dirname(config['ref_fasta']),
    threads: 4
    resources:
        time   = 2880,
        mem_mb = 12000
    shell:
        '''
            gautil coverage \
                {input.phase_vcf} \
                --refFolder={params.ref_dir}
        '''

