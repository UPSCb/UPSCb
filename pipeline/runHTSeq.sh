#!/bin/bash -l

## stop on error
set -ex

## usage
usage(){
echo >&2 \
"
	Usage: runHTSeq.sh [options] <out dir> <in.bam> <in.gff>	

	Options:
                -i precise the IDATTR
                   default to 'Parent', but e.g. should be 'pacid' 
                   for the P. trichocarpa gene exon gff3 file
                -s is the protocol stranded?
                   default to FALSE

        Note:
                BAM file are expected to be sorted by position
                Only HTSeq 0.6+ version(s) are supported
"
	exit 1
}

## Are we on UPPMAX?
if [ ! -z $SLURM_SUBMIT_DIR ]; then
    ## laod the modules
    echo Loading modules
    module load python/2.7.6
    module load bioinfo-tools
    module load samtools/0.1.19
else
    htseq=`which htseq-count`
    if [ "$?" -ne 0 ]; then
        echo "error: you need to install htseq or add it to your path"
        exit 1
    fi
fi

## check the version
isVersion6=`htseq-count --help | grep "version 0.6" |  wc -l`
if [ $isVersion6 != 1 ]; then
    echo Only HTSeq version 0.6+ are supported
    usage
fi

## options
IDATTR="Parent"
stranded=0

## get the options
while getopts i:s option
do
        case "$option" in
	    i) IDATTR=$OPTARG;;
	    s) stranded=1;;
	    \?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

## we get two dir and two files as input
if [ $# == 4 ]; then
    echo "This function arguments have changed!"
    usage
fi


if [ $# != 3 ]; then
    echo "This function takes one directory, one bam and one gff3 file as arguments"
    usage
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
fi

if [ ! -f $2 ]; then
    echo "The third argument needs to be an existing bam file"
    usage
fi
nam=`basename ${2//.bam/}`

if [ ! -f $3 ]; then
    echo "The forth argument needs to be an existing gff3 file"
    usage
fi

## sort by id
## samtools sort -n $3 $2/${nam}-byname

## get the count table
if [ $stranded == 0 ]; then
## since we are not using strand specific, go for the union
    htseq-count -f bam -r pos -m union -s no -t exon -i $IDATTR $2 $3 > $1/$nam.txt
else
    htseq-count -f bam -r pos -m intersection-nonempty -s reverse -t exon -i $IDATTR $2 $3 > $1/$nam.txt
fi

## clean
## rm $2/${nam}-byname.bam


