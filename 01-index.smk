import os

# Extract basename of fasta file (part between path and extension name)
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]

rule all:
    input:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa", "fai", "dict"]),
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"])

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
        mv {input.reference_genome}.fai {output}
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

rule hisat2_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.1.ht2",
        "genome_index/{ref_basename}.2.ht2",
        "genome_index/{ref_basename}.3.ht2",
        "genome_index/{ref_basename}.4.ht2",
        "genome_index/{ref_basename}.5.ht2",
        "genome_index/{ref_basename}.6.ht2",
        "genome_index/{ref_basename}.7.ht2",
        "genome_index/{ref_basename}.8.ht2"
    log:
        "logs/index/hisat2_index_{ref_basename}.log"
    shell:
        """
        hisat2-build \
        {input.reference_genome} \
        genome_index/{wildcards.ref_basename} \
        &> {log}
        """