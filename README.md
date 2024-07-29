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

### Input 

#### BAM Files
Example code for running SV-MeCa starting from BAM files and using reference genome build hg38. In this example BAM files are on the host machine under /inputbam, the reference file under /inputref, the bed file under /inputbed, and the user want to save the output under /output 

```
docker run -v /inputbam:/bam -v /inputref:/ref -v /inputbed:/bed -v /output:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh bam -bam /bam/yourbam.bam -ref /ref/yourref.fasta -sample yoursamplename -build hg38 -has_chr true" -bed /bed/yourbed.bed 
```
with `/workspace/SV-MeCa/results` static!

**Hint:** 
- Be aware of the correct use of quotation marks `"`.
- Be aware of the option `-has_chr`: If your data have chromosomes in the form `chr17`, it should be `true`, otherwise `false`
- Depending on the location of your files, fewer/more bind (-v) are possible. Docker doesn't allow more than one host path to be bound to a container path. 


#### VCF Files 

Example code for running SV-MeCa starting from VCF files and using reference genome build hg19. In this example VCF files and the statistic file are all into the same folder on the host machine under /your/directory/input/path and the user want to save the output under /your/directory/output/path

```
docker run -v /your/directory/input/path:/input -v /your/directory/output/path:/workspace/SV-MeCa/results wembasop/sv-meca:1.0 "/workspace/SV-MeCa/run_svmeca.sh vcf -bd /input/breakdancer.vcf -dl /input/delly.vcf -is /input/insurveyor.vcf -lp /input/lumpy.vcf -mt /input/manta.vcf -pd /input/pindel.vcf -td /input/tardis.vcf -st /input/stats.txt -sample yoursamplename -build hg19 -has_chr true" 
```
with `/workspace/SV-MeCa/results` static!

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
After the execution of SV-MeCa, a folder named by the provided sample name will be findable into the output folder.

The following structure should be expected:

Output_folder\
  - Sample name\
    - merge\ 
      - `samplename`.`sv_caller`.edit.sort.vcf : processed & filtered vcf files used for the merging step with SURVIVOR (Only `INS` and `DEL` )
      - `samplename`.merge.stats: some infos about the data: Amount of SVs per tool that have been merged by SURVIVOR, mean coverage, length of the reads
      - `samplename`.survivor.vcf: merged VCF file from SURVIVOR
    - metrics\
      - `samplename`.metrics: GATK CollectWgsMetrics output file 
    - standalone\
      - `samplename`.`sv_caller`.vcf: VCF files from the standalone SV callers (+ `.cfg` file from breakdancer)
      - `samplename`.breakdancer.ctx: Breakdancer original CTX output file 
    - sv-meca\
      - `samplename`.svmeca.vcf: SV-MeCa output containing `DEL` & `INS` SVs. The `QUAL` column contains scores in a range from [0-1000] corresponding to the probability (`QUAL = Prob *1000`) to represent true positive SVs.
    - report_`samplename`.html: [Nextflow report](https://www.nextflow.io/docs/latest/tracing.html#execution-report) containing metrics about the execution
    - timeline_`samplename`.html: [Nextflow timeline](https://www.nextflow.io/docs/latest/tracing.html#timeline-report) execution of each process
    - trace_`samplename`.tsv : [Nextflow trace](https://www.nextflow.io/docs/latest/tracing.html#trace-report) containing informations about each process executed in the pipeline: task_id,	hash,	native_id,	name,	status,	attempt,	exit,	realtime,	cpus,	%cpu,	memory,	%mem,	rss,	vmem,	peak_rss,	peak_vmem

#### SV-MeCa VCF

An example SV call in the VCF looks as follow:
```
1       991955  BD_2        C       <DEL>   306     LowProb IDS=BD_2,vh_del_2668;SUPP=2;SUPP_VEC=1000001;SVLEN=755;SVTYPE=DEL;END=992469
```

The SV call (Deletion) of length `755` on chromosome `1` at position `991955` with the ID `BD2`, has a score of `306`,i.e a probability of 0.306 to represent a TP SV. Since it's under 0.5 the filter is not PASS and the SV is flag with `LowProb`. The call is detected by two callers (`SUPP=2`), which are BreakDancer and TARDIS(`SUPP_VEC=1000001`). `SUPP_VEC` encodes in a binary form the standalone SV caller, which detected the putative SV, always in the following order: BreakDancer(1), Delly(0), INSurVeyor(0), Lumpy(0), Manta(0), Pindel(0), TARDIS(1). The INFO tag `IDS` contains the ID of the SV Call in the original output of the standalone SV callers([see](#folder-structure)). 

## SV-MeCa on HPC

### Singularity 
To run SV-MeCa with singularity you need a machine with root rigths. These steps should be done:

Machine with root rights
- Fetch the docker image as described in [installation](#installation)
- Build a singularity container from the docker-daemon as described [here](https://stackoverflow.com/questions/60314664/how-to-build-singularity-container-from-dockerfile)

HPC Machine / local machine (without root rights) 
- Copy the image to the HPC machine
- Run SV-MeCa with singularity considering the volumes to bind. 

### Clone GitHub repository
There is a way to execute SV-MeCa without docker which need advanced configuration requirements and enable the user to use different executors like SLURM or SGE. Here is list of [supported executors](https://www.nextflow.io/docs/latest/executor.html) 
Due to the amount of axecutor it's not possible to produce a comprehensive documentation. Feel free to open an issue or contact rudel.nkouamedjo-fankep@uk-koeln.de for further support on that question.

## Questions:

Feel free to open an issue for further questions, thanks!
