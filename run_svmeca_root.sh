#!/usr/bin/bash 
set -Eeuo pipefail

LOG_D=".log"
CONDA=""
SAMPLE=""
BAM=""
WORK="scratch/${SAMPLE}"
BED=""
REF=""
LOG="${LOG_D}/.nextflow_${SAMPLE}.log"
SONIC=""
RUN_OPT=""


export NXF_DISABLE_CHECK_LATEST=true

if [ ! -d "$LOG_D" ]; then
    # If the LOG_D doesn't exist, create it
    mkdir -p "$LOG_D"
fi

echo """
################ RUNNING COMMAND WITH PARAMETERS #############
nextflow -bg -log $LOG run main.nf \
        -w $WORK \
        --myConda $CONDA \
        --bed $BED \
        --input $BAM \
        --sample $SAMPLE \
        --sonic $SONIC \
        --containerBind "$WORK" \
        --reference $REF $RUN_OPT&>>"$LOG_D/${SAMPLE}_out_err.log"
################    END    ###################################

""" >"$LOG_D/${SAMPLE}_out_err.log"

nextflow -bg -log $LOG run main.nf \
        -w $WORK \
        --myConda $CONDA \
        --bed $BED \
        --input $BAM \
        --sample $SAMPLE \
        --sonic $SONIC \
        --containerBind "$WORK" \
        --reference $REF $RUN_OPT&>>"$LOG_D/${SAMPLE}_out_err.log"


unset NXF_DISABLE_CHECK_LATEST