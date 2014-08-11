#!/bin/bash -l

## stop on error but be verbose
set -e
set -x

##
#  Run fastQC 
##
## Usage:  sh runFastQC.sh file ouputFolder

## sanity checks
## executable
## are we on UPPMAX
if [ ! -z $SLURM_SUBMIT_DIR ]; then
	module load bioinfo-tools
	module load FastQC/0.10.1 
##	echo "Running on UPPMAX"
else
##	echo "Running locally"
	fastqc=`which fastqc`
	if [ "$?" == "1" ]; then
		echo "please install fastqc before running this script or add it to your PATH"
		exit 1
	fi

	if [ ! -f $fastqc -a ! -x $fastqc ]; then
		echo "your fastQC does not appear to be an executable file"
		exit 1
	fi
fi

## arguments
if [ $# != 2 ]; then
   echo "This script takes two arguments: the input file and the output directory"
   exit 1
fi

## input file
if [ ! -f $1 ]; then
	echo "The first argument needs to be an existing fastq (optionally gz) file"
	exit 1
fi

## output dir
if [ ! -d $2 ]; then
	echo "The second argument needs to be an existing output directory." 
fi

## start
fastqc --noextract --outdir $2 $1
