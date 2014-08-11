#!/bin/bash -l

usage() {
    echo "usage: `basename $0` <fastq>

Run fastQValidator on a FASTQ file. Prints output on stdout and
exits with a non-zero exit status if the input file does not
conform to the standard.

ARGUMENTS:
    fastq   a FASTQ file, can be gzipped

NOTES:
    fastQValidator must lie in your PATH" 1>&2
}

## stop on error
set -e

## check
if [ $# != 1 ]; then
    echo "The argument should be one fastq filename" 1>&2
    usage
    exit 1
fi

if [ ! -f $1 ]; then
    echo "The fastq filename you provided does not exist" 1>&2
    usage
    exit 1
fi

if ! hash fastQValidator 2>/dev/null; then
    echo "fastQValidator was not found in your path" 1>&2
    exit 1
fi

## we print 1000 errors, should be enough
fastQValidator --file $1 --printableErrors 1000
