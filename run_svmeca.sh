#!/bin/bash

# Script version
SCRIPT_VERSION="1.0"

LOG_D="/tmp/.log"

FMT=""
SAMPLE=""
BUILD=""
NO_CHR=""

INPUT=""
REF=""
BED=""

BD=""
DL=""
IS=""
LP=""
MT=""
PD=""
TD=""
STS=""
EXTRA=""

PARAMS=""

export NXF_DISABLE_CHECK_LATEST=true

if [ ! -d "$LOG_D" ]; then
    # If the LOG_D doesn't exist, create it
    mkdir -p "$LOG_D"
fi

usage() {
    echo "Usage: $0 <BAM|VCF> [options]"
    echo "Version: $SCRIPT_VERSION"
    echo "Description: This script is used to run the SV-MeCa pipeline,"
    echo "             meta-caller for structural variant detection," 
    echo "             depending on the user input format. BAM or VCF files."
    echo "Options:"
    echo "  For BAM:"
    echo "    -bam <bam>: Input BAM file"
    echo "    -ref <fasta|fa>: Reference genome file"
    echo "    -bed <bed>: BED file (optional)"
    echo "  For VCF:"
    echo "    -bd <vcf>: BreakDancer VCF"
    echo "    -dl <vcf>: Delly VCF"
    echo "    -is <vcf>: INSurVeyor VCF"
    echo "    -lp <vcf>: Lumpy VCF"
    echo "    -mt <vcf>: Manta VCF"
    echo "    -pd <vcf>: Pindel VCF"
    echo "    -td <vcf>: TARDIS VCF"
    echo "    -st <txt>: statistic txt file"
    echo "General options:"
    echo "    -sample <value>: Sample name (required)"
    echo "    -build <hg38|hg19>: Build value (required)"
    echo "    -no_chr <true|false>: Specify whether to include chr prefix (required)"
    echo "    -extra <extra_value>: Extra nextflow options (optional)"
    exit 1
}

if [ "$#" -lt 1 ]; then
    usage
fi

FMT=$1
shift

case $FMT in
    "bam")
        while [ "$#" -gt 0 ]; do
            case $1 in
                -bam)
                    INPUT=$2
                    shift 2
                    ;;
                -ref)
                    REF=$2
                    shift 2
                    ;;
                -bed)
                    BED=$2
                    shift 2
                    ;;
                -sample)
                    SAMPLE=$2
                    shift 2
                    ;;
                -build)
                    case $2 in
                        "hg38"|"hg19")
                            BUILD=$2
                            shift 2
                            ;;
                        *)
                            echo "Invalid value for -build. Accepted values are 'hg38' or 'hg19'"
                            usage
                            ;;
                    esac
                    ;;
                -no_chr)
                    case $2 in
                        "true"|"false")
                            NO_CHR=$2
                            shift 2
                            ;;
                        *)
                            echo "Invalid value for -no_chr. Accepted values are 'true' or 'false'"
                            usage
                            ;;
                    esac
                    ;;
                -extra)
                    EXTRA=" $2"
                    shift 2
                    ;;
                *)
                    usage
                    ;;
            esac
        done
        if [ -z "$INPUT" ] || [ -z "$REF" ] || [ -z "$SAMPLE" ] || [ -z "$BUILD" ] || [ -z "$NO_CHR" ]; then
            echo "Missing required parameters for BAM format"
            usage
        fi
        PARAMS="--input $INPUT --reference $REF --bed $BED" #"--build $BUILD --sample $SAMPLE --no_chr $NO_CHR"
        ;;
    "vcf")
        while [ "$#" -gt 0 ]; do
            case $1 in
                -bd)
                    BD=$2
                    shift 2
                    ;;
                -dl)
                    DL=$2
                    shift 2
                    ;;
                -is)
                    IS=$2
                    shift 2
                    ;;
                -lp)
                    LP=$2
                    shift 2
                    ;;
                -mt)
                    MT=$2
                    shift 2
                    ;;
                -pd)
                    PD=$2
                    shift 2
                    ;;
                -td)
                    TD=$2
                    shift 2
                    ;;
                -st)
                    STS=$2
                    shift 2
                    ;;
                -sample)
                    SAMPLE=$2
                    shift 2
                    ;;
                -build)
                    case $2 in
                        "hg38"|"hg19")
                            BUILD=$2
                            shift 2
                            ;;
                        *)
                            echo "Invalid value for -build. Accepted values are 'hg38' or 'hg19'"
                            usage
                            ;;
                    esac
                    ;;
                -no_chr)
                    case $2 in
                        "true"|"false")
                            NO_CHR=$2
                            shift 2
                            ;;
                        *)
                            echo "Invalid value for -no_chr. Accepted values are 'true' or 'false'"
                            usage
                            ;;
                    esac
                    ;;
                -extra)
                    EXTRA=" $2"
                    shift 2
                    ;;
                *)
                    usage
                    ;;
            esac
        done
        if [ -z "$BD" ] || [ -z "$DL" ] || [ -z "$IS" ] || [ -z "$LP" ] || [ -z "$MT" ] || [ -z "$PD" ] || [ -z "$TD" ] || [ -z "$STS" ] || [ -z "$SAMPLE" ] || [ -z "$BUILD" ] || [ -z "$NO_CHR" ]; then
            echo "Missing required parameters for VCF format"
            usage
        fi
        PARAMS="--breakdancer $BD --delly $DL --insurveyor $IS --lumpy $LP --manta $MT --pindel $PD --tardis $TD --stats $STS"
        ;;
    *)
        usage
        ;;
esac

LOG="${LOG_D}/.nextflow_${SAMPLE}.log"
WORK="/tmp/scratch/${SAMPLE}"
BASE="-w $WORK --containerBind $WORK --format $FMT"
PARAMS="--build $BUILD --sample $SAMPLE --no_chr $NO_CHR ${PARAMS}"

# Run nextflow command
echo ""
echo "nextflow run /workspace/SV-MeCa/main.nf $BASE $PARAMS$EXTRA"
echo ""
nextflow run /workspace/SV-MeCa/main.nf $BASE $PARAMS$EXTRA
