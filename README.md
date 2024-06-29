# SNP calling snakemake pipeline
> This pipeline is inspired by Yang Fanjing (Zhejiang University).

This pipeline is comprised by following part:
- Genome file index creation (BWA/samtools/GATK)
- Resequencing reads quality control
- Resequencing reads map to reference genome
- GATK SNP and INDEL calling pipeline
- SNP quality control (MissingRate, MAF etc.)

Main result file:
