import os

# Extract basename of fasta file (part between path and extension name)
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]

rule all:
    input:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa", "fai", "dict"])
        
rule bwa_index:
    input:
        reference_genome=config["ref"]
    output:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["amb", "ann", "bwt", "pac", "sa"])
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
