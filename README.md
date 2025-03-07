# SV-MeCa: Structural Variant Meta-Caller 


<!--ts-->
  * [Introduction](#introduction)
    * [Quality Metrics Considered for Deletion Calls (Summarizing Overview)](#quality-metrics-considered-for-deletion-calls-summarizing-overview)
    * [Quality Metrics Considered for Insertion & Duplication Calls (Summarizing Overview)](#quality-metrics-considered-for-insertion--duplication-calls-summarizing-overview)
  * [Implementation](#implementation)
  * [Installation](#installation)
  * [Basic Usage](#basic-usage)
    * [Help](#help)
    * [Input](#input)
      * [BAM File](#bam-file)
      * [VCF Files](#vcf-files)
    * [Output](#output)
      * [Folder Structure](#folder-structure)
      * [SV-MeCa VCF](#sv-meca-vcf)
    * [SV-MeCa on HPC](#sv-meca-on-hpc)
      * [Singularity](#singularity)
      * [Clone GitHub Repository](#clone-github-repository)
  * [Questions](#questions)
<!--te-->

## Introduction

SV-MeCa is a meta-caller for WGS short read data that combines seven standalone structural variant (SV) callers, including [BreakDancer](https://github.com/genome/breakdancer), [Delly](https://github.com/dellytools/delly), [INSurVeyor](https://github.com/kensung-lab/INSurVeyor), [LUMPY](https://github.com/arq5x/lumpy-sv), [Manta](https://github.com/Illumina/manta), [Pindel](https://github.com/genome/pindel), and [TARDIS](https://github.com/BilkentCompGen/tardis). The results from each caller are merged using the tool [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), considering the type of SV but not the strand.

What sets SV-MeCa apart from other meta-callers is its use of specific quality metrics from each SV caller (see Tables below). These metrics, normalized by coverage or read length if necessary, were used as inputs to train XGBoost decision tree classifiers that predict the probability of an SV call being a true positive, i.e., a real-wolrd variant.

### Quality Metrics Considered for Deletion Calls (Summarizing Overview)

| Quality Metric | BreakDancer | Delly | LUMPY | Manta | Pindel | TARDIS |
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

A minimal example workflow for testing your installation is available in our [SV-MeCA_data](https://github.com/ccfboc-bioinformatics/SV-MeCa_data/tree/main/test) repository.

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

### Input 

#### BAM File
Example code for running SV-MeCa starting from BAM files and using reference genome build hg38. In this example, BAM files are located under `/inputbam`, the reference FASTA file under `/inputref`, the optional input BED file under `/inputbed`, and the user wants to save the output under `/output`:

```
docker run -v /inputbam:/bam -v /inputref:/ref -v /inputbed:/bed -v /output:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh bam -bam /bam/yourbam.bam -ref /ref/yourref.fasta -sample yoursamplename -build hg38 -has_chr true" -bed /bed/yourbed.bed 
```
Note that `/workspace/SV-MeCa/results` is static.

**Hint:** 
- Be aware of the correct use of quotation marks `"`.
- Be aware of the option `-has_chr`: If your data have chromosomes in the form `chr17`, it should be `true`, otherwise `false`
- Depending on the location of your files, fewer/more bind (-v) are possible. Docker doesn't allow more than one host path to be bound to a container path. 


#### VCF Files 

Example code for running SV-MeCa starting from VCF files and using reference genome build hg19. In this example, VCF files and the statistic file are all located in the same folder named `/your/directory/input/path`, and the user wants to save the output under `/your/directory/output/path`:

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh vcf -bd /input/breakdancer.vcf -dl /input/delly.vcf -is /input/insurveyor.vcf -lp /input/lumpy.vcf -mt /input/manta.vcf -pd /input/pindel.vcf -td /input/tardis.vcf -st /input/stats.txt -sample yoursamplename -build hg19 -has_chr true" 
```
Note that `/workspace/SV-MeCa/results` is static.

**Hint:** 
- Be aware of the correct use of quotation marks `"`.
- Be aware of the option `-has_chr`: If your data have chromosomes in the form `chr17`, it should be `true`, otherwise `false`
- Depending on the location of your files, fewer/more bind (-v) are possible. Docker doesn't allow more than one host path to be bound to a container path. 
- All VCF files must be provided.

In VCF mode, the additional input file `stats.txt` is required to provide information about mean sequencing coverage and read lengths. 

[Example file](https://github.com/wembaSop/SV-MeCa/tree/master/Test/stats.txt) content:

```
coverage=50
readlen=65
```

Note that in BAM mode, information on sequencing coverage and read lengths are automatically fetched from input data. 

### Output

#### Folder Structure

SV-MeCa creates the following folder structure as output, with the parent directory named by `<sample_name>` as specified by the `-sample` argument: 

output_folder\
  - `<sample_name>\`
    - `merge\` 
      - `<sample_name>.sv_caller.edit.sort.vcf` : processed & filtered VCF files used for the merging step with SURVIVOR (contains exclusively insertion and deletion calls)
      - `<sample_name>.merge.stats` : some information about the input data: amount of SV calls per tool that have been merged by SURVIVOR, mean coverage, and length of the reads
      - `<sample_name>.survivor.vcf` : merged VCF file obtained from SURVIVOR
    - `metrics\`
      - `<samplename>.metrics` : GATK CollectWgsMetrics output file 
    - standalone\
      - `<sample_name>.sv_caller.vcf` : VCF files obtained from the standalone SV callers
      - `<sample_name>.breakdancer.ctx`, `<sample_name>.breakdancer.cfg` : raw BreakDancer output files 
    - `sv-meca\`
      - `<sample_name>.svmeca.vcf` : SV-MeCa output containing deletion (`DEL`) & insertion (`INS`) SV calls (see description below) 
    - `report_<sample_name>.html` : [Nextflow report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) containing metrics about the execution
    - `timeline_<sample_name.html` : [Nextflow timeline](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) execution of each process
    - `trace_<sample_name>.tsv` : [Nextflow trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) containing informations about each process executed in the pipeline: `task_id`,	`hash`,	`native_id`, `name`, `status`,	`attempt`, `exit`, `realtime`, `cpus`,	`%cpu`,	`memory`,	`%mem`,	`rss`, `vmem`,	`peak_rss`,	`peak_vmem`

#### SV-MeCa VCF

An example SV call in the SV-MeCa's output VCF `<sample_name>.svmeca.vcf` looks as follows:
```
1       991955  BD_2        C       <DEL>   306     LowProb IDS=BD_2,vh_del_2668;SUPP=2;SUPP_VEC=1000001;SVLEN=755;SVTYPE=DEL;END=992469
```

This is the call of a deletion (`SVTYPE`) of length 755bp (`SVLEN`), starting from position `991955` at chromosme `1`. 
The `QUAL` column (6th column) encodes evidence for true positive SV calls in a range from [0-1000]. Values refer to model-derived prediction probabilities multiplied by 1000. 
In the example shown, `QUAL` score `306` refers to a probability of 0.306 of this deletion to represent a true positive SV call, i.e., real-world SV. 
As the probability is below 0.5, the `LowProb` flag is set in the `FILTER` column.

The shown VCF entry has been consistently called by 2 SV callers (`SUPP=2`), namely BreakDancer and TARDIS, as `SUPP_VEC` encodes caller agreement as binary flags in the following order: 
1. BreakDancer
2.  Delly
3.  INSurVeyor
4.  LUMPY
5.  Manta
6.  Pindel
7.  TARDIS

The id `BD_2` in the `ID` column refers to the corresponding entry in the original SURVIVOR-generated VCF file in output folder  `<sample_name>\merge`, the ids given by `IDS` in the `INFO` colunm (here `BD_2` and `vh_del_2668`) refer to variant call ids in standalone caller VCFs in folder `<sample_name>\standalone\` (see section Folder Structure). 

## SV-MeCa on HPC

### Singularity 
To run SV-MeCa with Singularity you need a machine with root privileges. The following steps are required:

1. On a machine providing root privileges:
- Fetch the Docker image as described in [installation](#installation)
- Build a Singularity container from the Docker daemon as described [here](https://stackoverflow.com/questions/60314664/how-to-build-singularity-container-from-dockerfile).
2. on HPC or local machine (root privileges are not required): 
- Copy the image to the machine.
- Run SV-MeCa with Singularity considering the volumes to bind. 

### Clone GitHub Repository
There is a way to execute SV-MeCa without Docker which need advanced configuration requirements and enable the user to use different executors like SLURM or SGE. Here is list of [supported executors](https://www.nextflow.io/docs/latest/executor.html) 
Due to the amount of possible executors it's not feasible to provide a comprehensive documentation for all of them. Feel free to open an issue for further support on that question.

## Questions:

Feel free to open an issue in case of questions or remarks, thanks!

