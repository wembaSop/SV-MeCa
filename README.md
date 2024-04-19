# SV-MeCa: Structural Variant Meta Caller 


## Introduction

SV-MeCa is a meta-caller for WGS short read data that combines seven standalone structural variant (SV) callers, including [BreakDancer](https://github.com/genome/breakdancer), [Delly](https://github.com/dellytools/delly), [INSurVeyor](https://github.com/kensung-lab/INSurVeyor), [Lumpy](https://github.com/arq5x/lumpy-sv), [Manta](https://github.com/Illumina/manta), [Pindel](https://github.com/genome/pindel), and [TARDIS](https://github.com/BilkentCompGen/tardis). The results from each caller are merged using the tool [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), considering the type of SV but not the strand.

What sets SV-MeCa apart from other meta-callers is its use of specific quality metrics from each SV caller (see Tables below). These metrics, normalized by coverage or read length if necessary, were used as inputs to train XGBoost decision tree classifiers that predict the probability of an SV call being a true positive, i.e., a real-wolrd variant.

### Quality Metrics Considered for Deletion Calls (Summarizing Overview)

| Quality Metric | BreakDancer | Delly | Lumpy | Manta | Pindel | TARDIS |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
No. of Supporting Read Pairs | * | | |  * |  | * | 
No. of Supporting Split Reads | | * | | | | * |
Variant Fraction | | | * |  | * | |  |
Strand Bias | * |  | * | | | |
Quality assigned (e.g., QUAL value) | | * | | * |  | * |
Local Read Depth | | * | | | * | |
Split Read Consensus Alignment Quality | | * | | | | |  
Consensus Sequence Entropy | | * | | | | |
Microhomology Length | | * | | * | * |  |
Genotype Quality | | * | * | | | |
PRECISE vs IMPRECISE | | | * | | | |

### Quality Metrics Considered for Insertion & Duplication Calls (Summarizing Overview)

| Quality Metric | Delly | InSurVeyor | Manta | Pindel | TARDIS |
| --- | :---: | :---: | :---: | :---: | :---: | 
No. of Supporting Read Pairs |  | * | * |  |  |  
No. of Supporting Split Reads | | * | * | | * | 
Variant Fraction | | | |  | * |  |
Quality assigned (e.g., QUAL value) | * | | |  | * |
Local Read Depth | | * | | * | | 
Split Read Consensus Alignment Quality | * | | | | |
Consensus Sequence Entropy | * | | | | |
Microhomology Length | | * | * | * |   |
Genotype Quality | | | * | | | 
PRECISE vs IMPRECISE | | * | * | | |

## Implementation

SV-MeCa is implemented as a nextflow workflow and available as docker container for a better portability and user-friendly set up. 

## Installation

You only need the Docker engine installed to set up SV-MeCa. [see](https://docs.docker.com/engine/install/ubuntu/)

If Docker is already installed, SV-MeCa can be fetched from Docker Hub: https://hub.docker.com/r/wembasop/sv-meca 

## Basic Usage

### Help

The help of SV-MeCa give a gut overview of required files to run it.

Use this command to see the help:

`docker run wembasop/sv-meca:1.0 /workspace/SV-MeCa/run_svmeca.sh`

output:

```
Usage: ./run_svmeca.sh <bam|vcf> [options]
Version: 1.0
Description: This script is used to run the SV-MeCa pipeline,
             meta-caller for structural variant detection,
             depending on the user input format. BAM or VCF files.
Options:
  For bam:
    -bam <bam>: Input BAM file
    -ref <fasta|fa>: Reference genome file
    -bed <bed>: BED file (optional)
  For vcf:
    -bd <vcf>: BreakDancer VCF
    -dl <vcf>: Delly VCF
    -is <vcf>: INSurVeyor VCF
    -lp <vcf>: Lumpy VCF
    -mt <vcf>: Manta VCF
    -pd <vcf>: Pindel VCF
    -td <vcf>: TARDIS VCF
    -st <txt>: statistic txt file
General options:
    -sample <value>: Sample name (required)
    -build <hg38|hg19>: Build value (required)
    -no_chr <true|false>: Specify whether to include chr prefix (required)
    -extra <extra_value>: Extra nextflow options (optional)
```
SV-MeCa support currently two type of data, BAM files and VCF files. 

### For BAM files 
For the excution of BAM files the command could look like this for an execution with the genome build hg38:

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh bam -bam /input/yourbam.bam -ref /input/yourref.fasta -sample yoursamplename -build hg38 -no_chr true" 
```
* ' " ' are important to keep

### For VCF files 
For the excution of VCF files the command could look like this for an execution with the genome build hg19:

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh vcf -bd /input/breakdancer.vcf -dl /input/delly.vcf -is /input/insurveyor.vcf -lp /input/lumpy.vcf -mt /input/manta.vcf -pd /input/pindel.vcf -td /input/tardis.vcf -st /input/stats.txt -sample yoursamplename -build hg19 -no_chr true" 
```
* ' " ' are important to keep

The file `stats.txt` is required, as it provide informations abour coverage and readlen of the data. It looks like:

```
coverage=50
readlen=65
```

An example of it is provided [here](https://github.com/wembaSop/SV-MeCa/tree/master/Test/stats.txt) 

It's required only for the VCF mode, as the provided informations are automatically fetched from the data in the BAM modus.
