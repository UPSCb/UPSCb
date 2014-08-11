#!/bin/bash

# Preprocessing script for RNA-Seq data.

# THIS SCRIPT IS NOT TO BE RUN THROUGH SBATCH, USE BASH!

underline=`tput smul`
nounderline=`tput rmul`
bold=`tput bold`
normal=`tput sgr0`

usage() {
    echo "usage: bash `basename $0` [OPTIONS] <proj> <fastq1> <fastq2> <outdir>

Run the RNA-Seq preprocessing pipeline, i.e. FastQC, Trimmomatic, Sortmerna
and STAR. You throw in a pair of fastq files and it spits out a BAM file. Sweet!

${bold}ARGUMENTS:${normal}
    proj    project that this should be run as
    fastq1  forward fastq file
    fastq2  reverse fastq file
    outdir  directory to save everything in

${bold}OPTIONS:${normal}
    -h        show this help message and exit
    -S cfg    use SLURM for job queuing with given config file
    -T cfg    use TORQUE PBS for job queuing with given config file
    -s n      step at which to start (see ${underline}STEPS${nounderline})
    -e n      step at which to end (see ${underline}STEPS${nounderline})
    -g dir    path to STAR reference to use (required if STAR is included
              in pipeline)
    -G gff    gene model GFF3 file for STAR
    -H gff    gff3 file for ${underline}H${nounderline}TSeq, if not given, the star gff3 will be used
    -t        library is s${underline}t${nounderline}randed (currently only relevant for HTSeq)
    -i        IDATTR in GFF3 file to report counts (default: 'Parent')

${bold}STEPS:${normal}
    The steps of this script are as follows:

    1) fastQValidator
    2) FastQC
    3) SortMeRNA
    4) FastQC
    5) Trimmomatic
    6) FastQC
    7) STAR
    8) HTSeq-count

${bold}NOTES:${normal}
    * This script should not be run through sbatch, just use bash.
    * If the output directory already exists, content may be overwritten
    * If the SORTMERNA and UPSCb environmental variables don't exist, the
      script will guess them. Set them yourself just to be safe." 1>&2
}

# Check the number of arguments
if [ $# -lt 5 ]; then
    usage
    exit 1
fi

batch_system=
cfg_file=
pstart=1
pend=100
star_ref=
star_gff=
htseq_gff=
stranded=0
idattr="Parent"
# Parse the options
OPTIND=1
while getopts "hS:T:s:e:g:G:H:ti:" opt; do
    case "$opt" in
        h) usage; exit 1 ;;
        S) cfg_file=$OPTARG; batch_system="slurm" ;;
        T) cfg_file=$OPTARG; batch_system="torque" ;;
        s) pstart=$OPTARG ;;
        e) pend=$OPTARG ;;
        g) star_ref=$OPTARG ;;
        G) star_gff=$OPTARG ;;
        H) htseq_gff=$OPTARG ;;
        t) stranded=1 ;;
        i) idattr=$OPTARG ;;
        ?) usage; exit 1 ;;
    esac
done

shift $((OPTIND - 1))

# Check that valid batch system configuration was provided
if [ -z $batch_system ]; then
    echo "ERROR: Must select either SLURM or TORQUE PBS for job queuing, along with a configuration file" 1>&2
    exit 1
fi

if [ ! -f $cfg_file ]; then
    echo "ERROR: Configuration file for job queuing must exist - could not find file: '$cfg_file'" 1>&2
    exit 1
fi

# load configuration
source $cfg_file

# verify all configuration parameters have been loaded
cfg_params=( global_batch_args fastqc_batch_args fastqvalidator_batch_args 
    htseq_batch_args sortmerna_batch_args star_batch_args trimmomatic_batch_args )
    
for cfg_param in "${cfg_params[@]}"
do
    if [ -z "${!cfg_param}" ]; then
        echo "ERROR: undefined configuration parameter: '$cfg_param'" 1>&2
        exit 1
    fi
done

# Check the step parameters
[[ $pstart =~ ^[0-9]+$ ]] || {
    echo "ERROR: '$pstart' is not a valid start value" 1>&2
    exit 1
}
[[ $pend =~ ^[0-9]+$ ]] || {
    echo "ERROR: '$pend' is not a valid end value" 1>&2
    exit 1
}
[[ $pend -lt $pstart  ]] && {
    echo "ERROR: you can't finish before you start" 1>&2
    exit 1
}

# Check if STAR is included and then if the reference is set
if [ $pstart -le 7 ] && [ $pend -ge 7 ]; then
    if [ -z $star_ref ]; then
        echo >&2 "ERROR: You are running STAR but have not given a STAR reference"
        exit 1
    elif [ ! -d $star_ref ]; then
        echo >&2 "ERROR: Could not find STAR reference"
        exit 1
    fi
    if [ ! -z $star_gff ] && [ ! -f $star_gff ]; then
        echo >&2 "ERROR: Could not find gff file: '$star_gff'"
        exit 1
    fi
fi

# Check that HTSeq will get a GFF3 file (if it's included in the pipeline)
if [ $pstart -le 8 ] && [ $pend -ge 8 ]; then
    if [ -z $htseq_gff ] && [ -z $star_gff ]; then
        echo >&2 "ERROR: HTSeq needs a GFF3 file"
        exit 1
    elif [ -z $htseq_gff ]; then
        htseq_gff=$star_gff
    elif [ ! -z $htseq_gff ] && [ ! -f $htseq_gff ]; then
        echo >&2 "ERROR: Could not find gff file: '$htseq_gff'"
        exit 1
    fi
fi

# This variable holds the absolute path of this script, i.e. the
# repository's pipeline directory.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -e

proj=$1

if [ ! -f $2 ]; then
    echo "ERROR: fastq1 is not a file: '$2'" 1>&2
    usage
    exit 1
fi

if [ ! -f $3 ]; then
    echo "ERROR: fastq2 is not a file: '$3'" 1>&2
    usage
    exit 1
fi

fastq1=$2
fastq2=$3

# Sample name to use for output
s_prefix=${fastq1%_1.f*q.gz}
sname=`basename $s_prefix`

if [ ! -d `dirname $4` ]; then
    echo "ERROR: could not find parent directory for output directory" 1>&2
    usage
    exit 1
fi

## Set up the directory structure
outdir=$4
[[ ! -d $outdir ]] && mkdir $outdir

## Use the directory name as a job identifier
JOBNAME=`basename $outdir`

## FastQC
fastqc_raw="$outdir/fastqc/raw"
[[ ! -d $fastqc_raw ]] && mkdir -p $fastqc_raw
fastqc_sortmerna="$outdir/fastqc/sortmerna"
[[ ! -d $fastqc_sortmerna ]] && mkdir $fastqc_sortmerna
fastqc_trimmomatic="$outdir/fastqc/trimmomatic"
[[ ! -d $fastqc_trimmomatic ]] && mkdir $fastqc_trimmomatic

## SortmeRNA
sortmerna="$outdir/sortmerna"
[[ ! -d $sortmerna ]] && mkdir $sortmerna

## Trimmomatic
trimmomatic="$outdir/trimmomatic"
[[ ! -d $trimmomatic ]] && mkdir $trimmomatic

## STAR
star="$outdir/star"
[[ ! -d $star ]] && mkdir $star

## HTSeq
htseq="$outdir/htseq"
[[ ! -d $htseq ]] && mkdir $htseq

## Export some variables
[[ -z $UPSCb ]] && export UPSCb=$PIPELINE_DIR/..
## I will just assume that the sortmerna data is symlinked in the repo
[[ -z $SORTMERNADIR ]] && export SORTMERNADIR=$PIPELINE_DIR/../data/sortmerna
if [ ! -e $SORTMERNADIR ]; then
    echo "ERROR: could not find the sortmerna data in $SORTMERNADIR" 1>&2
    exit 1
fi

## Functions ##
run_sbatch_usage() {
    echo "usage: run_sbatch OPTIONS <batch script> [<script param> ...]

    Run a batch script and echo the job id

    OPTIONS
      -e  stderr (required)
      -o  stdout (required)
      -J  job name
      -d  dependency" 1>&2
}

cleanup_scancel() {
    #jobs=$@
    if [ $# -gt 0 ]; then
        echo "Canceling already started jobs: $@" 1>&2
        scancel $@
    fi
    exit 3
}

cleanup_qdel() {
    if [ $# -gt 0 ]; then
        echo "Canceling already started jobs: $@" 1>&2
        qdel $@
    fi
    exit 3
}

cleanup() {
    case "$batch_system" in
        slurm) echo `cleanup_scancel $@`;;
        torque) echo `cleanup_qdel $@`;;
        *) echo $"ERROR: batch system value has invalid value - $batch_system"  1>&2; exit 1;;
    esac
}

run_sbatch() {
    # Start a batch script and echo the job id
    if [ $# -lt 3 ]; then
        run_sbatch_usage
        exit 1
    fi

    OPTIND=1

    log_path=""
    out_path=""
    dependency=""
    batch_args=""

    while getopts "J:e:o:d:a:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            a) batch_args="${!OPTARG}" ;;
            ?) run_sbatch_usage; exit 1 ;;
        esac
    done

    shift $((OPTIND-1))

    script=$1
    shift

    if [ -z $jobname ]; then
        jobname="${sname}.RNAPreproc.${script}"
    fi

    # Check that the output file paths are given
    if [ -z $log_path ] || [ -z $out_path ]; then
        run_sbatch_usage
        exit 1
    fi
    # Check that the output file directories exist
    if [ ! -d `dirname $log_path` ] || [ ! -d `dirname $out_path` ]; then
        echo "ERROR: stderr and stdout paths must exist" 1>&2
        run_sbatch_usage
        exit 1
    fi

    sbatch_options=
    if [ ! -z $dependency ]; then
        sbatch_options="-d $dependency"
    fi

    sbatch_echo=`sbatch -A "$proj" \
                "$global_batch_args" \
                "$batch_args" \
                -J "$jobname" \
                -e "$log_path" \
                -o "$out_path" \
                $sbatch_options \
                $PIPELINE_DIR/$script $@`

    if [ $? -ne 0 ]; then
        echo "ERROR: submission failed" 1>&2
        cleanup ${JOBIDS[*]}
    fi

    echo ${sbatch_echo//[^0-9]/}
}

run_qsub_usage() {
    echo "usage: run_qsub OPTIONS <batch script> [<script param> ...]

    Run a batch script and echo the job id

    OPTIONS
      -e  stderr (required)
      -o  stdout (required)
      -J  job name
      -d  dependency" 1>&2
}

run_qsub() {
    # Start a batch script and echo the job id
    if [ $# -lt 3 ]; then
        run_qsub_usage
        exit 1
    fi

    OPTIND=1

    log_path=""
    out_path=""
    dependency=""
    batch_args=""

    while getopts "J:e:o:d:a:" opt; do
        case "$opt" in
            J) jobname=$OPTARG ;;
            e) log_path=$OPTARG ;;
            o) out_path=$OPTARG ;;
            d) dependency=$OPTARG ;;
            a) batch_args="${!OPTARG}" ;;
            ?) run_qsub_usage; exit 1 ;;
        esac
    done

    shift $((OPTIND-1))

    script=$1
    shift

    if [ -z $jobname ]; then
        jobname="${sname}.RNAPreproc.${script}"
    fi

    # Check that the output file paths are given
    if [ -z $log_path ] || [ -z $out_path ]; then
        run_qsub_usage
        exit 1
    fi
    # Check that the output file directories exist
    if [ ! -d `dirname $log_path` ] || [ ! -d `dirname $out_path` ]; then
        echo "ERROR: stderr and stdout paths must exist" 1>&2
        run_qsub_usage
        exit 1
    fi

    qsub_options=
    if [ ! -z $dependency ]; then
        qsub_options="-W depend=$dependency"
    fi

    params=$@
    qsub_echo=`echo "$PIPELINE_DIR/$script $params" | \
                qsub -A "$proj" \
                "$global_batch_args" \
                "$batch_args" \
                -N "$jobname" \
                -e "$log_path" \
                -o "$out_path" \
                $qsub_options \
                $tmp_script`

    if [ $? -ne 0 ]; then
        echo "ERROR: submission failed" 1>&2
        cleanup ${JOBIDS[*]}
    fi

    echo ${qsub_echo}
}

run_batch() {
    case "$batch_system" in
        slurm) echo `run_sbatch $@`;;
        torque) echo `run_qsub $@`;;
        *) echo $"ERROR: batch system value has invalid value - $batch_system" 1>&2; exit 1;;
    esac
}
## End functions ##

## Job ID array
JOBIDS=()

## Run fastQValidator
if [ $pstart -le 1 ]; then
    fastqv_id1=`run_batch \
        -a fastqvalidator_batch_args \
        -e $fastqc_raw/${sname}_1_validate.err \
        -o $fastqc_raw/${sname}_1_validate.out \
        -J ${sname}.RNAseq.FastQValidate1 \
        runFastQValidator.sh $fastq1`
    JOBIDS+=($fastqv_id1)

    fastqv_id2=`run_batch \
        -a fastqvalidator_batch_args \
        -e $fastqc_raw/${sname}_2_validate.err \
        -o $fastqc_raw/${sname}_2_validate.out \
        -J ${sname}.RNAseq.FastQCValidate2 \
        runFastQValidator.sh $fastq2`
    JOBIDS+=($fastqv_id2)
fi

## Run FastQC. Depends on a successful run of fastQValidator
if [ $pstart -le 2 ] && [ $pend -ge 2 ]; then
    dep1=
    dep2=
    if [ $pstart -lt 2 ]; then
        dep1="-d afterok:$fastqv_id1"
        dep2="-d afterok:$fastqv_id2"
    fi
    fqc_raw_id1=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_raw/${sname}_1_fastqc.err \
        -o $fastqc_raw/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.raw1 \
        $dep1 runFastQC.sh $fastq1 $fastqc_raw`
    JOBIDS+=($fqc_raw_id1)

    fqc_raw_id2=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_raw/${sname}_2_fastqc.err \
        -o $fastqc_raw/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.raw2 \
        $dep2 runFastQC.sh $fastq2 $fastqc_raw`
    JOBIDS+=($fqc_raw_id2)
fi

#TODO: Run fastqc stats

## Run Sortmerna. Depends on a successful run of FastQC
if [ $pstart -le 3 ] && [ $pend -ge 3 ]; then
    dep=
    if [ $pstart -lt 3 ]; then
        dep="-d afterok:$fqc_raw_id1:$fqc_raw_id2"
    fi
    sortmerna_id=`run_batch \
        -a sortmerna_batch_args \
        -e $sortmerna/${sname}_sortmerna.err \
        -o $sortmerna/${sname}_sortmerna.out \
        $dep \
        -J ${sname}.RNAseq.SortMeRNA \
        runSortmerna.sh $sortmerna $SNIC_TMP $fastq1 $fastq2`
    JOBIDS+=($sortmerna_id)
fi

fastq_sort_1="$sortmerna/${sname}_sortmerna_1.fq.gz"
fastq_sort_2="$sortmerna/${sname}_sortmerna_2.fq.gz"

# Run FastQC. Depends on a successful run of SortMeRNA
if [ $pstart -le 4 ] && [ $pend -ge 4 ]; then
    dep=
    if [ $pstart -lt 4 ]; then
        dep="-d afterok:$sortmerna_id"
    elif [ ! -f $fastq_sort_1 ] || [ ! -f $fastq_sort_2 ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    fqc_sort_id1=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_sortmerna/${sname}_1_fastqc.err \
        -o $fastqc_sortmerna/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.SortMeRNA1 \
        $dep runFastQC.sh $fastq_sort_1 $fastqc_sortmerna`
    JOBIDS+=($fqc_sort_id1)

    fqc_sort_id2=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_sortmerna/${sname}_2_fastqc.err \
        -o $fastqc_sortmerna/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.SortMeRNA2 \
        $dep runFastQC.sh $fastq_sort_2 $fastqc_sortmerna`
    JOBIDS+=($fqc_sort_id2)
fi

## Run trimmomatic. Depends on a successful run of SortMeRNA
if [ $pstart -le 5 ] && [ $pend -ge 5 ]; then
    dep=
    if [ $pstart -lt 5 ]; then
        # SortMeRNA has to finish, if it was started
        dep="-d afterok:$sortmerna_id"
    elif [ ! -f $fastq_sort_1 ] || [ ! -f $fastq_sort_2 ]; then
        # If we start here, make sure the files exist
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    trimmomatic_id=`run_batch \
        -a trimmomatic_batch_args \
        -e $trimmomatic/${sname}_trimmomatic.err \
        -o $trimmomatic/${sname}_trimmomatic.log \
        -J ${sname}.RNAseq.Trimmomatic \
        $dep \
        runTrimmomatic.sh $fastq_sort_1 $fastq_sort_2 $trimmomatic SLIDINGWINDOW:7:30 MINLEN:50`
    JOBIDS+=($trimmomatic_id)
fi

## Trimmed fastq file paths
fastq_trimmed_1="$trimmomatic/${sname}_sortmerna_trimmomatic_1.fq.gz"
fastq_trimmed_2="$trimmomatic/${sname}_sortmerna_trimmomatic_2.fq.gz"

# Run FastQC. Depends on a successful run of Trimmomatic
if [ $pstart -le 6 ] && [ $pend -ge 6 ]; then
    dep=
    if [ $pstart -lt 6 ]; then
        dep="-d afterok:$trimmomatic_id"
    elif [ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi
    fqc_trimmed_id1=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_trimmomatic/${sname}_1_fastqc.err \
        -o $fastqc_trimmomatic/${sname}_1_fastqc.out \
        -J ${sname}.RNAseq.FastQC.Trimmomatic1 \
        $dep runFastQC.sh $fastq_trimmed_1 $fastqc_trimmomatic`
    JOBIDS+=($fqc_trimmed_id1)

    fqc_trimmed_id2=`run_batch \
        -a fastqc_batch_args \
        -e $fastqc_trimmomatic/${sname}_2_fastqc.err \
        -o $fastqc_trimmomatic/${sname}_2_fastqc.out \
        -J ${sname}.RNAseq.FastQC.Trimmomatic2 \
        $dep runFastQC.sh $fastq_trimmed_2 $fastqc_trimmomatic`
    JOBIDS+=($fqc_trimmed_id2)
fi

# Run STAR. Depends on a successful run of Trimmomatic.
if [ $pstart -le 7 ] && [ $pend -ge 7 ]; then
    dep=
    if [ $pstart -lt 7 ]; then
        dep="-d afterok:$trimmomatic_id"
    elif [ ! -f $fastq_trimmed_1 ] || [ ! -f $fastq_trimmed_2 ]; then
        echo >&2 "ERROR: rRNA-filtered FASTQ-files could not be found"
        cleanup
    fi

    if [ ! -z $star_gff ]; then
        star_id=`run_batch \
            -a star_batch_args \
            -e $star/${sname}_STAR.err \
            -o $star/${sname}_STAR.out \
            -J ${sname}.RNAseq.STAR \
            $dep runSTAR.sh -o $star $fastq_trimmed_1 $fastq_trimmed_2 $star_ref $star_gff -- --outReadsUnmapped Fastx`
    else
        star_id=`run_batch \
            -a star_batch_args \
            -e $star/${sname}_STAR.err \
            -o $star/${sname}_STAR.out \
            -J ${sname}.RNAseq.STAR \
            $dep runSTAR.sh -o $star -g $fastq_trimmed_1 $fastq_trimmed_2 $star_ref -- --outReadsUnmapped Fastx`
    fi
    JOBIDS+=($star_id)
fi

# Run HTSeq. Depends on a successful run of STAR
if [ $pstart -le 8 ] && [ $pend -ge 8 ]; then
    # Load a proper python version
    module load python/2.7.6
    if hash htseq-count; then
        if ! htseq-count --help 2>&1 | grep "version 0.6" >/dev/null 2>&1; then
            echo >&2 "ERROR: HTSeq v0.6 or higher is required"
            cleanup
        fi
    else
        echo >&2 "ERROR: Could not find HTSeq in path"
        cleanup
    fi

    dep=
    if [ $pstart -lt 8 ]; then
        dep="-d afterok:$star_id"
    elif [ ! -f $star/${sname}_sortmerna_trimmomatic_STAR.bam ]; then
        echo >&2 "ERROR: STAR alignment BAM file could not be found"
        cleanup
    fi

    strand_arg=
    if [ $stranded -eq 1 ]; then
        strand_arg="-s"
    fi

    htseq_id=`run_batch \
        -a htseq_batch_args \
        -e $htseq/${sname}_HTSeq.err \
        -o $htseq/${sname}_HTSeq.out \
        -J ${sname}.RNAseq.HTSeq \
        $dep runHTSeq.sh -i $idattr $strand_arg $htseq $star/${sname}_sortmerna_trimmomatic_STAR.bam $htseq_gff`
    JOBIDS+=($htseq_id)
fi

echo "Successfully submitted ${#JOBIDS[@]} jobs for ${sname}: ${JOBIDS[@]}"

