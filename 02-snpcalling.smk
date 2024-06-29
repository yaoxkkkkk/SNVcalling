# Load config file and parameters
ref = config["ref"]
QD_snp = config["QD_snp"]
MQ_snp = config["MQ_snp"]
FS_snp = config["FS_snp"]
SOR_snp = config["SOR_snp"]
MQRS_snp = config["MQRS_snp"]
RPRS_snp = config["RPRS_snp"]
QUAL_snp = config["QUAL_snp"]
QD_indel = config["QD_indel"]
FS_indel = config["FS_indel"]
SOR_indel = config["SOR_indel"]
RPRS_indel = config["RPRS_indel"]
QUAL_indel = config["QUAL_indel"]
missingrate = config["missingrate"]
maf = config["maf"]

# Target file
rule all:
    input:
        expand("mapping/{sample}.sorted.markdup.bam",sample=config["sample"]),
        expand("vcf/gvcf/{sample}.g.vcf.gz",sample=config["sample"]),
        "vcf/raw.vcf.gz",
        "vcf/gvcf/vcf.list",
        "vcf/snp/filtered.snp.vcf.gz",
        "vcf/snp/clean.maf.snp",
        "vcf/indel/filtered.indel.vcf.gz",
        "vcf/filtered.vcf.gz",
        "vcf/clean"

# rule fastp:
    # input:
        # "raw_data/{sample}_1.fastq.gz",
        # "raw_data/{sample}_2.fastq.gz"
    # output:
        # "clean_data/{sample}.1_clean.fq.gz",
        # "clean_data/{sample}.2_clean.fq.gz",
        # "clean_data/{sample}.fastp.html"
    # log:
        # "logs/fastp/{sample}.log"
    # shell:
        # "fastp \
        # --thread 6 \
        # -i {input[0]} \
        # -I {input[1]} \
        # -o {output[0]} \
        # -O {output[1]} \
        # -h {output[2]} \
        # &> {log}"

# rule bwa_map:
    # input:
        # "clean_data/{sample}_trim_clean_1.fq.gz",
        # "clean_data/{sample}_trim_clean_2.fq.gz",
        # index_1=f"{ref}.amb",
        # index_2=f"{ref}.ann",
        # index_3=f"{ref}.bwt",
        # index_4=f"{ref}.pac",
        # index_5=f"{ref}.sa"
    # output:
        # temp("mapping/{sample}.sorted.bam")
    # threads:4
    # params:
        # rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    # log:
        # "logs/bwa/bwa_map_{sample}.log"
    # shell:
        # """
        # bwa mem \
        # -R '{params.rg}' \
        # -t {threads} \
        # {ref} {input[0]} {input[1]} \
        # | samtools view -Sb \
        # | samtools sort > {output} \
		# &> {log}
        # """

# rule RemoveDuplicates:
    # input:
        # input_bam="mapping/{sample}.sorted.bam"
    # output:
        # output_bam="mapping/{sample}.sorted.markdup.bam",
        # metrics_txt="mapping/{sample}.sorted.markdup_metrics.txt"
    # log:
        # "logs/RemoveDuplicates_{sample}.log"
    # shell:
        # """
        # gatk MarkDuplicates \
        # -I {input} \
        # -O {output[0]} \
        # -M {output[1]} \
        # --REMOVE_DUPLICATES true \
        # --CREATE_INDEX true \
        # &> {log}
        # """

rule HaplotypeCaller:
    input:
        "mapping/{sample}.sorted.markdup.bam"       
    output:
        "vcf/gvcf/{sample}.g.vcf.gz"
    log:
        "logs/vcf/gvcf/{sample}.gvcf.log"
    threads:4
    resources: mem_mb=32768
    shell:
        """
        gatk HaplotypeCaller \
        -R {ref} \
        -I {input[0]} \
        -ERC GVCF \
        -O {output} \
        &> {log}
        """

rule ExtractVCFlist:
    input:
        expand("vcf/gvcf/{sample}.g.vcf.gz",sample=config["sample"])
    output:
        "vcf/gvcf/vcf.list"
    params:
        gvcf_dir="vcf/gvcf/"
    shell:
        "find {params.gvcf_dir} -name '*.gz' > {output}"
        
rule ConsolidateGVCFs:
    input:
        "vcf/gvcf/vcf.list"
    output:
        "vcf/cohort.g.vcf.gz"
    log:
        "logs/vcf/ConsolidateGVCFs.log"
    shell:
        """
        gatk CombineGVCFs  \
        -R {ref} \
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
        "logs/vcf/GenotypeGVCFs.log"
    shell:
        """
        gatk GenotypeGVCFs \
        -R {ref} \
        -V {input} \
        -O {output} \
        &> {log}
        """

rule select_snp:
    input:
        "vcf/raw.vcf.gz"
    output:
        "vcf/snp/raw.snp.vcf.gz"
    log:
        "logs/vcf/snp.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -O {output} \
        --select-type-to-include SNP \
        &> {log}
        """

rule mark_snp:
    input:
        "vcf/snp/raw.snp.vcf.gz"
    output:
        "vcf/snp/filter.snp.vcf.gz"
    log:
        "logs/vcf/snp.filter.vcf.log"
    shell:
        """
        gatk VariantFiltration \
        -R {ref} \
        -V {input} \
        --filter-expression "QD < {QD_snp}" \
        --filter-name 'SNP_QD_filter' \
        --filter-expression "MQ < {MQ_snp}" \
        --filter-name 'SNP_MQ_filter' \
        --filter-expression "FS > {FS_snp}" \
        --filter-name 'SNP_FS_filter' \
        --filter-expression "SOR > {SOR_snp}" \
        --filter-name 'SNP_SOR_filter' \
        --filter-expression "MQRankSum < {MQRS_snp}" \
        --filter-name 'SNP_MQRS_filter' \
        --filter-expression "ReadPosRankSum < {RPRS_snp}" \
        --filter-name 'SNP_RPRS_filter' \
        --filter-expression "QUAL < {QUAL_snp}" \
        --filter-name 'SNP_QUAL_filter' \
        -O {output} \
         &> {log}
        """
        
rule filter_snp:
    input:
        "vcf/snp/filter.snp.vcf.gz"
    output:
        "vcf/snp/filtered.snp.vcf.gz"
    log:
        "logs/vcf/snp.filtered.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -R {ref} \
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
        "logs/vcf/clean.maf.snp.vcf.log"
    shell:
        """
        vcftools \
        --gzvcf {input} \
        --max-missing {missingrate} \
        --maf {maf} \
        --out {output} \
        --recode \
         &> {log}
        """
        
rule select_indel:
    input:
        "vcf/raw.vcf.gz"
    output:
        "vcf/indel/raw.indel.vcf.gz"
    log:
        "logs/vcf/indel.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -O {output} \
        --select-type-to-include INDEL \
        &> {log}
        """
        
rule mark_indel:
    input:
        "vcf/indel/raw.indel.vcf.gz"
    output:
        "vcf/indel/filter.indel.vcf.gz"
    log:
        "logs/vcf/indel.filter.vcf.log"
    shell:
        """
        gatk VariantFiltration \
        -R {ref} \
        -V {input} \
        --filter-expression "QD < {QD_indel}" \
        --filter-name 'INDEL_QD_filter' \
        --filter-expression "FS > {FS_indel}" \
        --filter-name 'INDEL_FS_filter' \
        --filter-expression "SOR > {SOR_indel}" \
        --filter-name 'INDEL_SOR_filter' \
        --filter-expression "ReadPosRankSum < {RPRS_indel}" \
        --filter-name 'INDEL_RPRS_filter' \
        --filter-expression "QUAL < {QUAL_indel}" \
        --filter-name 'INDEL_QUAL_filter' \
        -O {output} \
         &> {log}
        """
        
rule filter_indel:
    input:
        "vcf/indel/filter.indel.vcf.gz"
    output:
        "vcf/indel/filtered.indel.vcf.gz"
    log:
        "logs/vcf/indel.filtered.vcf.log"
    shell:
        """
        gatk SelectVariants \
        -R {ref} \
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
        "logs/vcf/filtered.vcf.log"
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
        "logs/vcf/clean.vcf.log"
    shell:
        """
        vcftools \
        --gzvcf {input} \
        --max-missing {missingrate} \
        --out {output} \
        --recode \
         &> {log}
        """
