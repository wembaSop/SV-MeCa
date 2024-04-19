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

## Basic Usage
