# SV-MeCa: Structural Variant Meta-Caller 


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

| Quality Metric | Delly | INSurVeyor | Manta | Pindel | TARDIS |
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

Installation of the [Docker Engine](https://docs.docker.com/engine/install/) is required to set up SV-MeCa.

If Docker is installed, SV-MeCa can be fetched from Docker Hub [here](https://hub.docker.com/r/wembasop/sv-meca). 

## Basic Usage

### Help

The help of SV-MeCa (available via `docker run wembasop/sv-meca:1.0 /workspace/SV-MeCa/run_svmeca.sh`)  provides an overview of required input files:

```
Usage: ./run_svmeca.sh <bam|vcf> [options]
Version: 1.0
Description: This script is used to run the SV-MeCa pipeline,
             meta-caller for structural variant detection,
             depending on the user input format. BAM or VCF files.
Options:
  For bam input:
    -bam <bam>: Input BAM file
    -ref <fasta|fa>: Reference genome file
    -bed <bed>: BED file (optional)
  For vcf input:
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
    -has_chr <true|false>: Specify whether the input chromosome data has chr prefix (required)
    -extra <extra_value>: Extra nextflow options (optional)
```
SV-MeCa supports two types of initial input data: BAM files and VCF files. 

### BAM Files 

Example code for running SV-MeCa starting from BAM files and using reference genome build hg38:

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh bam -bam /input/yourbam.bam -ref /input/yourref.fasta -sample yoursamplename -build hg38 -has_chr true" 
```

**Hint:** 
- Be aware of the correct use of quotation marks `"`.
- Be aware of the option `-has_chr`: If your data have chromosomes in the form `chr17`, it should be `true`, otherwise `false`


### VCF Files 

Example code for running SV-MeCa starting from VCF files and using reference genome build hg19:

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh vcf -bd /input/breakdancer.vcf -dl /input/delly.vcf -is /input/insurveyor.vcf -lp /input/lumpy.vcf -mt /input/manta.vcf -pd /input/pindel.vcf -td /input/tardis.vcf -st /input/stats.txt -sample yoursamplename -build hg19 -has_chr true" 
```

**Hint:** 
- Be aware of the correct use of quotation marks `"`.
- All VCF files must be provided.

In VCF mode, the additional input file `stats.txt` is required to provide information about mean sequencing coverage and read lengths. 

[Example file](https://github.com/wembaSop/SV-MeCa/tree/master/Test/stats.txt) content:

```
coverage=50
readlen=65
```

Note that in BAM mode, information on sequencing coverage and read lengths are automatically fetched from input data. 
