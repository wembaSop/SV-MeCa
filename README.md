# SV-MeCa: Structural Variant Meta Caller 


## Introduction

SV-MeCa is a meta-caller that combines seven standalone Structural Variant (SV) callers, including Breakdancer, Delly, Insurveyor, Lumpy, Manta, Pindel, and TARDIS. The results from each caller are merged using the tool SURVIVOR, considering the type of SV but not the strand.

What sets SV-MeCa apart from other meta-callers is its use of specific quality metrics from each SV caller (see Table X). These metrics, normalized by coverage or read length if necessary, are used as inputs to train a machine learning model that predicts the probability of an SV call being a true positive event

### Quality Metrics Considered for Deletion Calls (Summarizing Overview)

| Quality Metric | BreakDancer | Delly | Lumpy | Manta | Pindel | TARDIS |
| --- | :---: | :---: | :---: | :---: | :---: | :---: |
No. of Supporting Read Pairs | * | | |  * |  | * | 
No. of Supproting Split Reads | | * | | | | * |
Variant Fraction | | | * |  | * | |  |
Strand Bias | * |  | * | | | |
Quality assigned (e.g., QUAL value) | | * | | * |  | * |
Local Read Depth | | * | | | * | |
Split Read Consensus Alignment Quality | | * | | | | |  
Consensus Sequence Entropy | | * | | | | |
Microhomology Length | | * | | * | * |  |
Genotype Quality | | * | * | | | |
PRECISE vs IMPRECISE | | | * | | | |

## Implementation

SV-MeCa is implemented as a nextflow workflow and available as docker container for a better portability and user-friendly set up. 

## Installation

## Basic Usage
