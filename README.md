# Population joint SNP calling snakemake pipeline

## Dependent Software

- fastp
- BWA
- GATK
- samtools
- VCFtools

## What the pipeline does

- Genome file index creation (BWA/samtools/GATK)
- Resequencing reads quality control
- Resequencing reads map to reference genome
- GATK SNP and INDEL calling pipeline
- SNP quality control (Missing rate, MAF etc.)

## What to input

- Species reference genome
- Population resequence data

## What to output

## Usage

### 1. Prepare your working directory

```shell
├── script
│   └── snake_pipeline
├── raw_data
├── genome_index
└── logs
```

Please storage your resequence data in `raw_data/` folder and genome file in `genome_index/` folder. Script files, pipeline files and configuration files can be stored in the way you are used to.

### 2. Prepare the config file

2.1 Move the genome file to `genome_index/` folder and add the genome fasta file absolute path like:

```shell
# Absolute path to the genome fasta file
ref: "/workingdir/genome_index/genome.fasta" 
```

2.2 Sometimes the fastq files may be ended with `.fastq.gz` or `.fq.gz`, specify the suffix of the fastq files if it's necessary.

```shell
# Fastq file suffix
fastq_suffix: " " # Default value is ".fq.gz"
```

2.3 Fill in the name of the samples. The samples name need to be filled with specific format like:

```shell
# Sample list
sample:
    - "1_5_001"
    - "1_5_002"
    - "1_5_003"
    - "1_5_004"
    - "1_5_005"
    - "1_5_006"
```

### 3. Submit the pipeline to HPC cluster

For example:

```shell
snakemake \
	--snakefile ${working_dir}/00-script/snake_pipeline/${snakemake_file} \
	-d ${working_dir} \
    --config ${working_dir}/SNPcalling_config.yaml \
	--cores ${cores_num} \
	--rerun-incomplete \
	--latency-wait 360 \
	--keep-going
```