import os

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

qualified_quality_phred = config.get("qualified_quality_phred", 20)
unqualified_percent_limit = config.get("unqualified_percent_limit", 40)
trim_front = config.get("trim_front", 10)

rule all:
    input:
        expand("mapping/{sample}.sorted.markdup.bam", sample=config["sample"]),
        expand("vcf/gvcf/{sample}.g.vcf.gz", sample=config["sample"]),
        expand("vcf/gvcf/{sample}.g.vcf.gz.tbi", sample=config["sample"]),
        "vcf/raw.vcf.gz",
        "vcf/gvcf/vcf.list",
        "vcf/snp/filtered.snp.vcf.gz",
        "vcf/snp/clean.maf.snp",
        "vcf/indel/filtered.indel.vcf.gz",
        "vcf/filtered.vcf.gz",
        "vcf/clean"

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
        r2="clean_data/{sample}_2_clean.fq.gz"
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
        "logs/gatk/RemoveDuplicates_{sample}.log"
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

rule HaplotypeCaller:
    input:
        reference_genome=config["ref"],
        bam="mapping/{sample}.sorted.markdup.bam"
    output:
        "vcf/gvcf/{sample}.g.vcf.gz",
        "vcf/gvcf/{sample}.g.vcf.gz.tbi"
    log:
        "logs/gatk/vcf/gvcf/{sample}.gvcf.log"
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
        -O {output} \
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

rule ConsolidateGVCFs:
    input:
        "vcf/gvcf/vcf.list"
    output:
        "vcf/cohort.g.vcf.gz"
    log:
        "logs/gatk/vcf/ConsolidateGVCFs.log"
    shell:
        """
        gatk CombineGVCFs  \
        -R {config['ref']} \
        -V {input} \
        -O {output} \
        &> {log}
        """

rule GenotypeGVCFs:
    input:
        "vcf/cohort.g.vcf.gz"
    output:
        "vcf/raw.vcf.gz"
    log:
        "logs/gatk/vcf/GenotypeGVCFs.log"
    shell:
        """
        gatk GenotypeGVCFs \
        -R {config['ref']} \
        -V {input} \
        -O {output} \
        &> {log}
        """

rule Select_snp:
    input:
        "vcf/raw.vcf.gz"
    output:
        "vcf/snp/raw.snp.vcf.gz"
    log:
        "logs/gatk/vcf/snp.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -O {output} \
        --select-type-to-include SNP \
        &> {log}
        """

rule Mark_snp:
    input:
        "vcf/snp/raw.snp.vcf.gz"
    output:
        "vcf/snp/filter.snp.vcf.gz"
    log:
        "logs/gatk/vcf/snp.filter.vcf.log"
    shell:
        """
        gatk VariantFiltration \
        -R {config['ref']} \
        -V {input} \
        --filter-expression "QD < {config['QD_snp']}" \
        --filter-name 'SNP_QD_filter' \
        --filter-expression "MQ < {config['MQ_snp']}" \
        --filter-name 'SNP_MQ_filter' \
        --filter-expression "FS > {config['FS_snp']}" \
        --filter-name 'SNP_FS_filter' \
        --filter-expression "SOR > {config['SOR_snp']}" \
        --filter-name 'SNP_SOR_filter' \
        --filter-expression "MQRankSum < {config['MQRS_snp']}" \
        --filter-name 'SNP_MQRS_filter' \
        --filter-expression "ReadPosRankSum < {config['RPRS_snp']}" \
        --filter-name 'SNP_RPRS_filter' \
        --filter-expression "QUAL < {config['QUAL_snp']}" \
        --filter-name 'SNP_QUAL_filter' \
        -O {output} \
        &> {log}
        """

rule Filter_snp:
    input:
        "vcf/snp/filter.snp.vcf.gz"
    output:
        "vcf/snp/filtered.snp.vcf.gz"
    log:
        "logs/gatk/vcf/snp.filtered.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -R {config['ref']} \
        -V {input} \
        --exclude-filtered \
        -O {output} \
        &> {log}
        """

rule SNPMissingRateAndMAFFilter:
    input:
        "vcf/snp/filtered.snp.vcf.gz"
    output:
        "vcf/snp/clean.maf.snp"
    log:
        "logs/gatk/vcf/clean.maf.snp.vcf.log"
    shell:
        """
        vcftools \
        --gzvcf {input} \
        --max-missing {config['missingrate']} \
        --maf {config['maf']} \
        --out {output} \
        --recode \
        &> {log}
        """

rule Select_indel:
    input:
        "vcf/raw.vcf.gz"
    output:
        "vcf/indel/raw.indel.vcf.gz"
    log:
        "logs/gatk/vcf/indel.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -O {output} \
        --select-type-to-include INDEL \
        &> {log}
        """

rule Mark_indel:
    input:
        "vcf/indel/raw.indel.vcf.gz"
    output:
        "vcf/indel/filter.indel.vcf.gz"
    log:
        "logs/gatk/vcf/indel.filter.vcf.log"
    shell:
        """
        gatk VariantFiltration \
        -R {config['ref']} \
        -V {input} \
        --filter-expression "QD < {config['QD_indel']}" \
        --filter-name 'INDEL_QD_filter' \
        --filter-expression "FS > {config['FS_indel']}" \
        --filter-name 'INDEL_FS_filter' \
        --filter-expression "SOR > {config['SOR_indel']}" \
        --filter-name 'INDEL_SOR_filter' \
        --filter-expression "ReadPosRankSum < {config['RPRS_indel']}" \
        --filter-name 'INDEL_RPRS_filter' \
        --filter-expression "QUAL < {config['QUAL_indel']}" \
        --filter-name 'INDEL_QUAL_filter' \
        -O {output} \
        &> {log}
        """

rule Filter_indel:
    input:
        "vcf/indel/filter.indel.vcf.gz"
    output:
        "vcf/indel/filtered.indel.vcf.gz"
    log:
        "logs/gatk/vcf/indel.filtered.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -R {config['ref']} \
        -V {input} \
        --exclude-filtered \
        -O {output} \
        &> {log}
        """

rule MergeSNPandINDEL:
    input:
        "vcf/snp/filtered.snp.vcf.gz",
        "vcf/indel/filtered.indel.vcf.gz"
    output:
        "vcf/filtered.vcf.gz"
    log:
        "logs/gatk/vcf/filtered.vcf.log"
    shell:
        """
        gatk MergeVcfs  \
        -I {input[0]} \
        -I {input[1]} \
        -O {output} \
        &> {log}
        """

rule VCFMissingRateFilter:
    input:
        "vcf/filtered.vcf.gz"
    output:
        "vcf/clean"
    log:
        "logs/gatk/vcf/clean.vcf.log"
    shell:
        """
        vcftools \
        --gzvcf {input} \
        --max-missing {config['missingrate']} \
        --out {output} \
        --recode \
        &> {log}
        """
