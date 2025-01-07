> This pipeline is inspired by F.J. Yang (Bioinplant Lab, Zhejiang University).

## Dependent Software

- [fastp](https://github.com/OpenGene/fastp)
- [BWA](https://github.com/lh3/bwa)
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [samtools](https://github.com/samtools/samtools)
- [VCFtools](https://github.com/vcftools/vcftools)
- bgzip (part of [htslib](https://github.com/samtools/htslib))
- [PanDepth](https://github.com/HuiyangYu/PanDepth)

## What to input

- Reference genome file
- WGS fastq files

## What to output

- SNV Basic set (missingrate < 0.6 and MAF > 0.005) and core set (missingrate < 0.8 and MAF > 0.05)
- SNP core set (no missing and MAF > 0.05, could be used for [population structure analysis](https://github.com/yaoxkkkkk/Population-structure-analysis-pipeline))

## Usage

### 1. Prepare your working directory

```shell
├── raw_data
├── genome_index
└── logs
```

Please storage your resequence data in `raw_data/` folder and genome file in `genome_index/` folder. Script files, pipeline files and configuration files can be stored the way you like.

### 2. Prepare the config file

The config file needs to be at the same folder of snakefile.

#### 2.1 Move the genome file to `genome_index/` folder and add the genome fasta file absolute path like:

```shell
# Absolute path to the genome fasta file
ref: "/workingdir/genome_index/genome.fasta" 
```

And chromosome IDs to genotype variants by per chromosomes to parallise the process.

```shell
chromosomes:
    - Chr01_hap1
    - Chr02_hap1
    - Chr03_hap1
    - ...
    - Chrnn_hap1
```

You can use the following command to generate the list.

```shell
cat /workingdir/genome_index/genome.fasta | grep ">" | awk '{gsub(/^>/, "    - "); print}'
```

#### 2.2 Sometimes the fastq files may be ended with `.fastq.gz` or `.fq.gz`, specify the suffix of the fastq files if it's necessary.

```shell
# Fastq file suffix
fastq_suffix: " " # Default value is ".fq.gz"
```

#### 2.3 Fill in the name of the samples. The samples name need to be filled with specific format like:

```shell
# Sample list, samples' name should start with letters.
sample:
    - "sample1"
    - "sample2"
    - "sample3"
    - "sample4"
    - ...
    - "samplen"
```

You can use following command to add sample list to the config file if you have a sample list txt file (for example `sample.list`):

```shell
# sample.list
sample1
sample2
sample3
sample4

# Add samples to the config file:
awk '{print "    - \"" $0 "\""}' sample.list >> ${working_dir}/SNPcalling_config.yaml
```

### 3. Submit the pipeline to HPC cluster

Put snakefile and configuration file in the same directory, then running it.

For example:

```shell
snakemake \
	--snakefile ${snakefile} \
	-d ${working_dir} \
	--cores ${cores_num}
```
