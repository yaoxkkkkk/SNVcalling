import os
import gzip

configfile: "SNPcalling_config.yaml"

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

qualified_quality_phred = config.get("qualified_quality_phred", 20)
unqualified_percent_limit = config.get("unqualified_percent_limit", 40)
trim_front = config.get("trim_front", 10)

rule all:
    input:
        "vcf/clean.vcf.gz",
        "vcf/snp/clean.maf.snp.vcf.gz",
        "mapping/merged_depth_stats.txt"

rule bwa_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.amb",
        "genome_index/{ref_basename}.ann",
        "genome_index/{ref_basename}.bwt",
        "genome_index/{ref_basename}.pac",
        "genome_index/{ref_basename}.sa"
    log:
        "logs/index/bwa_index_{ref_basename}.log"
    shell:
        """
        bwa index \
        -p genome_index/{wildcards.ref_basename} \
        {input.reference_genome} \
        &> {log}
        """
rule samtools_fai_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.fai"
    log:
        "logs/index/samtools_index_{ref_basename}.log"
    shell:
        """
        samtools faidx {input.reference_genome} &> {log}
        cp {input.reference_genome}.fai {output}
        """
rule gatk_dict_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.dict"
    log:
        "logs/index/gatk_index_{ref_basename}.log"
    shell:
        """
        gatk CreateSequenceDictionary \
        -R {input.reference_genome} \
        -O {output} \
        &> {log}
        """

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_1{fastq_suffix}",
        f"raw_data/{{sample}}_2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "logs/fastp/fastp_report/{sample}.fastp.html"
    threads: 2
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
        --thread {threads} \
        -i {input[0]} \
        -I {input[1]} \
        -o {output[0]} \
        -O {output[1]} \
        -h {output[2]} \
        -j /dev/null \
        -q {qualified_quality_phred} \
        -u {unqualified_percent_limit} \
        -f {trim_front} \
        &> {log}
        """

rule BWA_map:
    input:
        r1="clean_data/{sample}_1_clean.fq.gz",
        r2="clean_data/{sample}_2_clean.fq.gz",
        bwa_index=expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp("mapping/{sample}.sorted.bam")
    threads: 8
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA",
        prefix=f"genome_index/{ref_basename}"
    log:
        "logs/bwa/bwa_map_{sample}.log"
    shell:
        """
        bwa mem \
        -R '{params.rg}' \
        -t {threads} \
        {params.prefix} {input.r1} {input.r2} \
        | samtools sort -@ {threads} -o {output} \
        &> {log}
        """

rule RemoveDuplicates:
    input:
        "mapping/{sample}.sorted.bam"
    output:
        "mapping/{sample}.sorted.markdup.bam",
        "mapping/markdup_metrics/{sample}.sorted.markdup_metrics.txt"
    log:
        "logs/RemoveDuplicates_{sample}.log"
    shell:
        """
        gatk MarkDuplicates \
        -I {input} \
        -O {output[0]} \
        -M {output[1]} \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true \
        &> {log}
        """

rule BAMDepthStat:
    input:
        bam="mapping/{sample}.sorted.markdup.bam"
    output:
        temp("mapping/{sample}.chr.stat.gz")
    log:
        "logs/depth/{sample}_pandepth.log"
    threads: 2
    shell:
        """
        pandepth \
        -i {input.bam} \
        -o mapping/{wildcards.sample} \
        -t {threads} \
        &> {log}
        """

rule MergeDepthStats:
    input:
        expand("mapping/{sample}.chr.stat.gz", sample=config["sample"])
    output:
        "mapping/merged_depth_stats.txt"
    log:
        "logs/depth/merge_depth_stats.log"
    run:
        with open(output[0], 'w') as out_file:
            for stat_file in input:
                sample = stat_file.split("/")[-1].split(".")[0]
                with gzip.open(stat_file, 'rt') as f: 
                    last_line = f.readlines()[-1].strip()
                    out_file.write(f"{sample}\t{last_line}\n") 

rule HaplotypeCaller:
    input:
        reference_genome=config["ref"],
        bam="mapping/{sample}.sorted.markdup.bam",
        samtools_index="genome_index/{ref_basename}.fai",
        gatk_dict_index="genome_index/{ref_basename}.dict"
    output:
        "vcf/gvcf/{sample}.g.vcf.gz",
        "vcf/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        "logs/vcf/gvcf/{sample}.gvcf.log"
    threads: 4
    params:
        ERC="GVCF"
    resources:
        mem_mb=32768
    shell:
        """
        gatk HaplotypeCaller \
        -R {input.reference_genome} \
        -I {input.bam} \
        -ERC {params.ERC} \
        -O {output[0]} \
        &> {log}
        """

rule ExtractVCFlist:
    input:
        expand("vcf/gvcf/{sample}.g.vcf.gz", sample=config["sample"])
    output:
        "vcf/gvcf/vcf.list"
    params:
        gvcf_dir="vcf/gvcf/"
    shell:
        """
        find {params.gvcf_dir} -name '*.g.vcf.gz' > {output}
        """

rule ConsolidateGVCFsPerChromosome:
    input:
        reference_genome=config["ref"],
        vcflist="vcf/gvcf/vcf.list"
    output:
        temp("vcf/cohort_{chrom}.g.vcf.gz"),
        temp("vcf/cohort_{chrom}.g.vcf.gz.tbi")
    log:
        "logs/vcf/ConsolidateGVCFs_{chrom}.log"
    threads: 4
    shell:
        """
        gatk CombineGVCFs \
        -R {input.reference_genome} \
        -V {input.vcflist} \
        -L {wildcards.chrom} \
        -O {output[0]} \
        &> {log}
        """

rule GenotypeGVCFsPerChromosome:
    input:
        reference_genome=config["ref"],
        cohort_gvcf="vcf/cohort_{chrom}.g.vcf.gz"
    output:
        temp("vcf/raw_{chrom}.vcf.gz"),
        temp("vcf/raw_{chrom}.vcf.gz.tbi")
    log:
        "logs/vcf/GenotypeGVCFs_{chrom}.log"
    threads: 4
    resources:
        mem_mb=32768
    shell:
        """
        gatk GenotypeGVCFs \
        -R {input.reference_genome} \
        -V {input.cohort_gvcf} \
        -O {output[0]} \
        &> {log}
        """

rule SelectSNPsPerChromosome:
    input:
        raw_vcf="vcf/raw_{chrom}.vcf.gz"
    output:
        temp("vcf/snp/raw_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/raw_{chrom}.snp.vcf.gz.tbi")
    log:
        "logs/vcf/snp_{chrom}.log"
    threads: 4
    shell:
        """
        gatk SelectVariants \
        -V {input.raw_vcf} \
        -O {output[0]} \
        --select-type-to-include SNP \
        &> {log}
        """

rule MarkFilteredSNPsPerChromosome:
    input:
        reference_genome=config["ref"],
        snp_vcf="vcf/snp/raw_{chrom}.snp.vcf.gz"
    output:
        temp("vcf/snp/mark_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/mark_{chrom}.snp.vcf.gz.tbi")
    log:
        "logs/vcf/snp_mark_{chrom}.log"
    params:
        QD_snp=config["QD_snp"],
        MQ_snp=config["MQ_snp"],
        FS_snp=config["FS_snp"],
        SOR_snp=config["SOR_snp"],
        MQRS_snp=config["MQRS_snp"],
        RPRS_snp=config["RPRS_snp"],
        QUAL_snp=config["QUAL_snp"]
    shell:
        """
        gatk VariantFiltration \
        -R {input.reference_genome} \
        -V {input.snp_vcf} \
        --filter-expression 'QD < {params.QD_snp}' \
        --filter-name 'SNP_QD_filter' \
        --filter-expression 'MQ < {params.MQ_snp}' \
        --filter-name 'SNP_MQ_filter' \
        --filter-expression 'FS > {params.FS_snp}' \
        --filter-name 'SNP_FS_filter' \
        --filter-expression 'SOR > {params.SOR_snp}' \
        --filter-name 'SNP_SOR_filter' \
        --filter-expression 'MQRankSum < {params.MQRS_snp}' \
        --filter-name 'SNP_MQRS_filter' \
        --filter-expression 'ReadPosRankSum < {params.RPRS_snp}' \
        --filter-name 'SNP_RPRS_filter' \
        --filter-expression 'QUAL < {params.QUAL_snp}' \
        --filter-name 'SNP_QUAL_filter' \
        -O {output[0]} \
        &> {log}
        """

rule FilterSNPsPerChromosome:
    input:
        reference_genome=config["ref"],
        mark_snp_vcf="vcf/snp/mark_{chrom}.snp.vcf.gz"
    output:
        temp("vcf/snp/filter_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/filter_{chrom}.snp.vcf.gz.tbi")
    log:
        "logs/vcf/snp_filter_{chrom}.log"
    shell:
        """
        gatk SelectVariants \
        -R {input.reference_genome} \
        -V {input.mark_snp_vcf} \
        --exclude-filtered \
        -O {output[0]} \
        &> {log}
        """
		
rule ExtractSNPList:
    input:
        expand("vcf/snp/filter_{chrom}.snp.vcf.gz", chrom=config["chromosomes"])
    output:
        temp("vcf/snp/snp_vcf.list")
    params:
        vcf_dir="vcf/snp/"
    shell:
        """
        find {params.vcf_dir} -name 'filter_*.snp.vcf.gz' > {output}
        """
		
rule MergeSNPs:
    input:
        "vcf/snp/snp_vcf.list"
    output:
        "vcf/snp/filter.snp.vcf.gz",
        "vcf/snp/filter.snp.vcf.gz.tbi"
    log:
        "logs/vcf/merge_filter_snp.log"
    shell:
        """
        gatk MergeVcfs \
        -I {input} \
        -O {output[0]} \
        &> {log}
        """
		
rule SelectIndelsPerChromosome:
    input:
        raw_vcf=temp("vcf/raw_{chrom}.vcf.gz")
    output:
        temp("vcf/indel/raw_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/raw_{chrom}.indel.vcf.gz.tbi")
    log:
        "logs/vcf/indel_{chrom}.log"
    threads: 4
    shell:
        """
        gatk SelectVariants \
        -V {input.raw_vcf} \
        -O {output[0]} \
        --select-type-to-include INDEL \
        &> {log}
        """
		
rule MarkFilteredIndelsPerChromosome:
    input:
        reference_genome=config["ref"],
        indel_vcf="vcf/indel/raw_{chrom}.indel.vcf.gz"
    output:
        temp("vcf/indel/mark_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/mark_{chrom}.indel.vcf.gz.tbi")
    log:
        "logs/vcf/indel_mark_{chrom}.log"
    params:
        QD_indel=config["QD_indel"],
        FS_indel=config["FS_indel"],
        SOR_indel=config["SOR_indel"],
        RPRS_indel=config["RPRS_indel"],
        QUAL_indel=config["QUAL_indel"]
    shell:
        """
        gatk VariantFiltration \
        -R {input.reference_genome} \
        -V {input.indel_vcf} \
        --filter-expression 'QD < {params.QD_indel}' \
        --filter-name 'INDEL_QD_filter' \
        --filter-expression 'FS > {params.FS_indel}' \
        --filter-name 'INDEL_FS_filter' \
        --filter-expression 'SOR > {params.SOR_indel}' \
        --filter-name 'INDEL_SOR_filter' \
        --filter-expression 'ReadPosRankSum < {params.RPRS_indel}' \
        --filter-name 'INDEL_RPRS_filter' \
        --filter-expression 'QUAL < {params.QUAL_indel}' \
        --filter-name 'INDEL_QUAL_filter' \
        -O {output[0]} \
        &> {log}
        """

rule FilterIndelsPerChromosome:
    input:
        reference_genome=config["ref"],
        mark_indel_vcf="vcf/indel/mark_{chrom}.indel.vcf.gz"
    output:
        temp("vcf/indel/filter_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/filter_{chrom}.indel.vcf.gz.tbi")
    log:
        "logs/vcf/indel_clean_{chrom}.log"
    shell:
        """
        gatk SelectVariants \
        -R {input.reference_genome} \
        -V {input.mark_indel_vcf} \
        --exclude-filtered \
        -O {output[0]} \
        &> {log}
        """
		
rule ExtractINDELList:
    input:
        expand("vcf/indel/filter_{chrom}.indel.vcf.gz", chrom=config["chromosomes"])
    output:
        temp("vcf/indel/indel_vcf.list")
    params:
        vcf_dir="vcf/indel/"
    shell:
        """
        find {params.vcf_dir} -name 'filter_*.indel.vcf.gz' > {output}
        """
		
rule MergeINDELs:
    input:
        "vcf/indel/indel_vcf.list"
    output:
        "vcf/indel/filter.indel.vcf.gz",
        "vcf/indel/filter.indel.vcf.gz.tbi"
    log:
        "logs/vcf/merge_filter_indel.log"
    shell:
        """
        gatk MergeVcfs \
        -I {input} \
        -O {output[0]} \
        &> {log}
        """
		
rule MergeSNPandINDEL:
    input:
        snp_vcf="vcf/snp/filter.snp.vcf.gz",
        indel_vcf="vcf/indel/filter.indel.vcf.gz"
    output:
        "vcf/filter.vcf.gz",
        "vcf/filter.vcf.gz.tbi"
    log:
        "logs/vcf/filter.vcf.log"
    shell:
        """
        gatk MergeVcfs  \
        -I {input.snp_vcf} \
        -I {input.indel_vcf} \
        -O {output[0]} \
        &> {log}
        """

rule SNPMissingRateAndMAFFilter:
    input:
        filtered_snp_vcf="vcf/snp/filter.snp.vcf.gz"
    output:
        "vcf/snp/clean.maf.snp.vcf.gz"
    log:
        "logs/vcf/clean.maf.snp.vcf.log"
    params:
        missingrate=config["missingrate"],
        maf=config["maf"]
    threads: 8
    shell:
        """
        vcftools \
        --gzvcf {input.filtered_snp_vcf} \
        --max-missing {params.missingrate} \
        --maf {params.maf} \
        --recode \
        --stdout \
        | bgzip > {output} \
        &> {log}
        """
		
rule VCFMissingRateFilter:
    input:
        filtered_vcf="vcf/filter.vcf.gz"
    output:
        "vcf/clean.vcf.gz"
    log:
        "logs/vcf/clean.vcf.log"
    params:
        missingrate=config["missingrate"]
    threads: 8
    shell:
        """
        vcftools \
        --gzvcf {input.filtered_vcf} \
        --max-missing {params.missingrate} \
        --maf 0.0005 \
        --recode \
        --stdout \
        | bgzip > {output} \
        &> {log}
        """
