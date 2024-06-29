import os

# Load the config file
configfile: "SNPcalling_config.yaml"

# Extract basename of fasta file (part between path and extension name)
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]

rule all:
    input:
        expand("genome_index/{ref_basename}.amb", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.ann", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.bwt", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.pac", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.sa", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.fai", ref_basename=ref_basename),
        expand("genome_index/{ref_basename}.dict", ref_basename=ref_basename)

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

rule fai_index:
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

rule dict_index:
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
