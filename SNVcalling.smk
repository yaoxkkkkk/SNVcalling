import os
import gzip

configfile: "SNVcalling_config.yaml"

ref_basename=os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix=config.get("fastq_suffix")

rule all:
    input:
        expand("vcf/gvcf/{sample}.g.vcf.gz", sample=config["sample"])
        "mapping/merged_depth_stats.txt"
        "vcf/snv.core.vcf.gz"

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
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        bwa index \
        -p genome_index/{wildcards.ref_basename} \
        {input.reference_genome} \
        2> {log}
        """

rule samtools_fai_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.fai"
    log:
        "logs/index/samtools_index_{ref_basename}.log"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=2000
    shell:
        """
        samtools faidx {input.reference_genome} 2> {log}
        cp {input.reference_genome}.fai {output}
        """

rule gatk_dict_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.dict"
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/index/gatk_index_{ref_basename}.log"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=4000
    shell:
        """
        gatk CreateSequenceDictionary \
        -R {input.reference_genome} \
        -O {output} \
        2> {log}
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
    resources:
        cpus_per_task=2,
        mem_mb=4000
    params:
        qualified_quality_phred=config["qualified_quality_phred"],
        unqualified_percent_limit=config["unqualified_percent_limit"],
        trim_front=config["trim_front"]
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
        -q {params.qualified_quality_phred} \
        -u {params.unqualified_percent_limit} \
        -f {params.trim_front} \
        2> {log}
        """

rule BWA_map:
    input:
        r1="clean_data/{sample}_1_clean.fq.gz",
        r2="clean_data/{sample}_2_clean.fq.gz",
        bwa_index=expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp("mapping/{sample}.sorted.bam")
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=32000
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
        2> {log}
        """

rule RemoveDuplicates:
    input:
        "mapping/{sample}.sorted.bam"
    output:
        "mapping/{sample}.sorted.markdup.bam",
        "mapping/markdup_metrics/{sample}.sorted.markdup_metrics.txt"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/bwa/RemoveDuplicates_{sample}.log"
    shell:
        """
        gatk MarkDuplicates \
        -I {input} \
        -O {output[0]} \
        -M {output[1]} \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true \
        2> {log}
        """

rule BAMDepthStat:
    input:
        bam="mapping/{sample}.sorted.markdup.bam"
    output:
        temp("mapping/{sample}.chr.stat.gz")
    log:
        "logs/depth/{sample}_pandepth.log"
    threads: 2
    resources:
        cpus_per_task=2,
        mem_mb=4000
    shell:
        """
        pandepth \
        -i {input.bam} \
        -o mapping/{wildcards.sample} \
        -t {threads} \
        2> {log}
        """

rule MergeDepthStats:
    input:
        expand("mapping/{sample}.chr.stat.gz", sample=config["sample"])
    output:
        "mapping/merged_depth_stats.txt"
    log:
        "logs/depth/merge_depth_stats.log"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=2000
    run:
        with open(output[0], 'w') as out_file:
            for stat_file in input:
                sample=stat_file.split("/")[-1].split(".")[0]
                with gzip.open(stat_file, 'rt') as f: 
                    last_line=f.readlines()[-1].strip()
                    out_file.write(f"{sample}\t{last_line}\n")

rule HaplotypeCaller:
    input:
        reference_genome=config["ref"],
        bam="mapping/{sample}.sorted.markdup.bam",
        samtools_index=f"genome_index/{ref_basename}.fai",
        gatk_dict_index=f"genome_index/{ref_basename}.dict"
    output:
        "vcf/gvcf/{sample}.g.vcf.gz",
        "vcf/gvcf/{sample}.g.vcf.gz.tbi"
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/gvcf/{sample}.gvcf.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=16000
    params:
        ERC="GVCF"
    shell:
        """
        gatk HaplotypeCaller \
        -R {input.reference_genome} \
        -I {input.bam} \
        -ERC {params.ERC} \
        -O {output[0]} \
        2> {log}
        """

rule ExtractVCFlist:
    input:
        expand("vcf/gvcf/{sample}.g.vcf.gz", sample=config["sample"])
    output:
        "vcf/gvcf/vcf.list"
    params:
        gvcf_dir="vcf/gvcf/"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=1000
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
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/ConsolidateGVCFs_{chrom}.log"
    threads: 20
    resources:
        cpus_per_task=20,
        mem_mb=80000
    shell:
        """
        gatk CombineGVCFs \
        -R {input.reference_genome} \
        -V {input.vcflist} \
        -L {wildcards.chrom} \
        -O {output[0]} \
        2> {log}
        """

rule GenotypeGVCFsPerChromosome:
    input:
        reference_genome=config["ref"],
        cohort_gvcf="vcf/cohort_{chrom}.g.vcf.gz",
        cohort_gvcf_index="vcf/cohort_{chrom}.g.vcf.gz.tbi"
    output:
        temp("vcf/raw_{chrom}.vcf.gz"),
        temp("vcf/raw_{chrom}.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/GenotypeGVCFs_{chrom}.log"
    threads: 16
    resources:
        cpus_per_task=16,
        mem_mb=64000
    shell:
        """
        gatk GenotypeGVCFs \
        -R {input.reference_genome} \
        -V {input.cohort_gvcf} \
        -O {output[0]} \
        2> {log}
        """

rule SelectSNPsPerChromosome:
    input:
        raw_vcf="vcf/raw_{chrom}.vcf.gz",
        raw_vcf_index="vcf/raw_{chrom}.vcf.gz.tbi"
    output:
        temp("vcf/snp/raw_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/raw_{chrom}.snp.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/snp_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        gatk SelectVariants \
        -V {input.raw_vcf} \
        -O {output[0]} \
        --select-type-to-include SNP \
        2> {log}
        """

rule MarkFilteredSNPsPerChromosome:
    input:
        reference_genome=config["ref"],
        snp_vcf="vcf/snp/raw_{chrom}.snp.vcf.gz",
        snp_vcf_index="vcf/snp/raw_{chrom}.snp.vcf.gz.tbi"
    output:
        temp("vcf/snp/mark_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/mark_{chrom}.snp.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/snp_mark_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=12000
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
        2> {log}
        """

rule FilterSNPsPerChromosome:
    input:
        reference_genome=config["ref"],
        mark_snp_vcf="vcf/snp/mark_{chrom}.snp.vcf.gz",
        mark_snp_vcf_index="vcf/snp/mark_{chrom}.snp.vcf.gz.tbi"
    output:
        temp("vcf/snp/filter_{chrom}.snp.vcf.gz"),
        temp("vcf/snp/filter_{chrom}.snp.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/snp_filter_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        gatk SelectVariants \
        -R {input.reference_genome} \
        -V {input.mark_snp_vcf} \
        --exclude-filtered \
        -O {output[0]} \
        2> {log}
        """

rule SNPCoresetPerChromosome:
    input:
        "vcf/snp/filter_{chrom}.snp.vcf.gz",
        "vcf/snp/filter_{chrom}.snp.vcf.gz.tbi"
    output:
        temp("vcf/snp/coreset_{chrom}.snp.vcf.gz")
    log:
        "logs/vcf/snp_coreset_{chrom}.log"
    params:
        missingrate=config["core"]["missingrate"],
        maf=config["core"]["maf"]
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=16000
    shell:
        """
        vcftools \
        --gzvcf {input[0]} \
        --max-alleles 2 \
        --max-missing {params.missingrate} \
        --maf {params.maf} \
        --recode \
        --stdout \
        | bgzip > {output} \
        2> {log}
        """
        
rule SelectIndelsPerChromosome:
    input:
        raw_vcf="vcf/raw_{chrom}.vcf.gz",
        raw_vcf_index="vcf/raw_{chrom}.vcf.gz.tbi"
    output:
        temp("vcf/indel/raw_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/raw_{chrom}.indel.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/indel_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        gatk SelectVariants \
        -V {input.raw_vcf} \
        -O {output[0]} \
        --select-type-to-include INDEL \
        --max-indel-size 50 \
        2> {log}
        """
        
rule MarkFilteredIndelsPerChromosome:
    input:
        reference_genome=config["ref"],
        indel_vcf="vcf/indel/raw_{chrom}.indel.vcf.gz",
        indel_vcf_index="vcf/indel/raw_{chrom}.indel.vcf.gz.tbi"
    output:
        temp("vcf/indel/mark_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/mark_{chrom}.indel.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/indel_mark_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=12000
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
        2> {log}
        """

rule FilterIndelsPerChromosome:
    input:
        reference_genome=config["ref"],
        mark_indel_vcf="vcf/indel/mark_{chrom}.indel.vcf.gz",
        mark_indel_vcf_index="vcf/indel/mark_{chrom}.indel.vcf.gz.tbi"
    output:
        temp("vcf/indel/filter_{chrom}.indel.vcf.gz"),
        temp("vcf/indel/filter_{chrom}.indel.vcf.gz.tbi")
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/indel_clean_{chrom}.log"
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        gatk SelectVariants \
        -R {input.reference_genome} \
        -V {input.mark_indel_vcf} \
        --exclude-filtered \
        -O {output[0]} \
        2> {log}
        """

rule IndelCoresetPerChromosome:
    input:
        "vcf/indel/filter_{chrom}.indel.vcf.gz",
        "vcf/indel/filter_{chrom}.indel.vcf.gz.tbi"
    output:
        temp("vcf/indel/coreset_{chrom}.indel.vcf.gz")
    log:
        "logs/vcf/indel_coreset_{chrom}.log"
    params:
        missingrate=config["core"]["missingrate"],
        maf=config["core"]["maf"]
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=16000
    shell:
        """
        vcftools \
        --gzvcf {input[0]} \
        --max-alleles 2 \
        --max-missing {params.missingrate} \
        --maf {params.maf} \
        --recode \
        --stdout \
        | bgzip > {output} \
        2> {log}
        """

rule ExtractcoresetList:
    input:
        expand("vcf/snp/coreset_{chrom}.snp.vcf.gz", chrom=config["chromosomes"]),
        expand("vcf/indel/coreset_{chrom}.indel.vcf.gz", chrom=config["chromosomes"])
    output:
        temp("vcf/coreset_vcf.list")
    params:
        snp_vcf_dir="vcf/snp/",
        indel_vcf_dir="vcf/indel/"
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=1000
    shell:
        """
        find {params.snp_vcf_dir} -name 'coreset_*.snp.vcf.gz' > {output}
        find {params.indel_vcf_dir} -name 'coreset_*.indel.vcf.gz' >> {output}
        """

rule MergeCoreSet:
    input:
        "vcf/coreset_vcf.list",
        expand("vcf/snp/coreset_{chrom}.snp.vcf.gz", chrom=config["chromosomes"]),
        expand("vcf/indel/coreset_{chrom}.indel.vcf.gz", chrom=config["chromosomes"])
    output:
        "vcf/snv.core.vcf.gz",
        "vcf/snv.core.vcf.gz.tbi"
    conda:
        config["conda_env"]["gatk_env"]
    log:
        "logs/vcf/merge_coreset.log"
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=32000
    shell:
        """
        gatk MergeVcfs \
        -I {input[0]} \
        -O {output[0]} \
        2> {log}
        """
