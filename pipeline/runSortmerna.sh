#!/bin/bash -l
#SBATCH -p node
## for large files
## we don't need the proc but the mem
## we could give that as param
#SBATCH -n 16
## time too for large files
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
## mail-user and A have to be set in the submit script

## stop on error
set -e

## be verbose and extend the commands
set -x

## check the options if any
KEEP=0
useMtSSU=1
UNPAIRED=0
PROC=16

## local run
## replaced by checking for the SORTMERNADIR - see below
## if [ -z $SLURM_SUBMIT_DIR ]; then
##    SLURM_SUBMIT_DIR=`pwd`
## fi

## usage
usage(){
echo >&2 \
"
	Usage: runSortmerna.sh [option] <out dir> <tmp dir> <forward fastq.gz> <reverse fastq.gz>
	
	Options:
                -k keep the rRNA
                -m do not run against mtSSU
		-p number of threads to be used (default $PROC)
                -u single end data (in that case only the forward fastq is needed)

         Note:
               1) The SORTMERNADIR environment variable needs to be set
               2) Only SortMeRna version 1.9 is supported
"
	exit 1
}

## get the options
while getopts kmp:u option
do
        case "$option" in
	    k) KEEP=1;;
	    m) useMtSSU=0;;
	    p) PROC=$OPTARG;;
	    u) UNPAIRED=1;;
		\?) ## unknown flag
		usage;;
        esac
done
shift `expr $OPTIND - 1`

##
echo Setting up

## set some env var
## this location is not in Git anymore!
## it has to be downloaded by the user
## check the ethylene-insensitive project submitter to see
## how to set that up
if [ -z $SORTMERNADIR ]; then
    echo You need to set your SORTMERNADIR environment variable
    usage
fi

## set the dbs
db5s=$SORTMERNADIR/rRNA_databases/rfam-5s-database-id98.fasta
db58s=$SORTMERNADIR/rRNA_databases/rfam-5.8s-database-id98.fasta
db16s=$SORTMERNADIR/rRNA_databases/silva-bac-16s-database-id85.fasta
db18s=$SORTMERNADIR/rRNA_databases/silva-euk-18s-database-id95.fasta
db23s=$SORTMERNADIR/rRNA_databases/silva-bac-23s-database-id98.fasta
db28s=$SORTMERNADIR/rRNA_databases/silva-euk-28s-database-id98.fasta
dbNum=6
dbs="$db5s $db58s $db16s $db18s $db23s $db28s"
if [ $useMtSSU == 1 ]; then
    mtSSU=$SORTMERNADIR/rRNA_databases/mtSSU_UCLUST-95-identity.fasta
    dbs="$db5s $db58s $db16s $db18s $db23s $db28s $mtSSU"
    dbNum=7
fi



##
echo Checking

## we get two dir and two files as input
if [ $UNPAIRED == 0 ]; then
    if [ $# != 4 ]; then
	echo "This function takes two directories and two files as arguments"
	usage
    fi
else
    if [ $# != 3 ]; then
	echo "This function takes two directories and one file as argument"
	usage
    fi
fi

if [ ! -d $1 ]; then
    echo "The first argument needs to be an existing directory"
    usage
fi

if [ ! -d $2 ]; then
    echo "The second argument needs to be an existing directory"
    usage
fi

## 
echo Gunzipping

## unzip the files
if [ ! -f $3 ]; then
    echo "The third argument needs to be an existing fastq.gz file"
    usage
fi
f1=`basename ${3//.gz/}`

if [ $UNPAIRED == 0 ]; then
    if [ ! -f $4 ]; then
	echo "The forth argument needs to be an existing fastq.gz file"
	usage
    fi
    f2=`basename ${4//.gz/}`
fi

## decompress them
if [ ! -f $2/$f1 ]; then
    gunzip -c $3 > $2/$f1
fi
if [ $UNPAIRED == 0 ]; then
    if [ ! -f $2/$f2 ]; then
	gunzip -c $4 > $2/$f2
    fi
fi

## interleave them
fm=`basename ${3//.f*q.gz/}`
if [ $UNPAIRED == 0 ]; then
    isVersion9=`sortmerna --version | grep "version 1.9" | wc -l`
    if [ $isVersion9 != 1 ]; then
	echo Only SortMeRna version 1.9 is supported
	usage
    else
	merge-paired-reads.sh $2/$f1 $2/$f2 $2/$fm
    fi
fi

##
if [ $UNPAIRED == 0 ]; then
    echo Pre-cleaning
    rm -f $2/$f1 $2/$f2
fi

##
echo Sorting

## PE
if [ $UNPAIRED == 0 ]; then
    fo=`basename ${3//_[1,2].f*q.gz/_sortmerna}`
else
    fo=`basename ${3//.f*q.gz/_sortmerna}`
fi

## check the options
opt=
if [ $KEEP -eq 1 ]; then
    opt="--bydbs --accept $2/${fo}_rRNA"
fi 

if [ $UNPAIRED == 0 ]; then
    sortmerna -n $dbNum --db $dbs --I $2/$fm --other $2/$fo --log $1/$fo -a $PROC -v --paired-in $opt
else
    sortmerna -n $dbNum --db $dbs --I $2/$f1 --other $1/$fo --log $1/$fo -a $PROC -v $opt
fi

## deinterleave it
if [ $UNPAIRED == 0 ]; then
    ## sortmerna get confused by dots in the filenames
    if [ ! -f $2/$fo.fastq ]; then
	mv $2/$fo.* $2/$fo.fastq
    fi
    unmerge-paired-reads.sh $2/$fo.fastq $1/${fo}_1.fq $1/${fo}_2.fq
fi

## rm the tmp
echo Post-Cleaning
if [ $UNPAIRED == 0 ]; then
    rm -f $2/$fm $2/$fo.fastq
else
    rm -f $2/$f1
fi

## deinterleave the rest if needed
if [ $KEEP -eq 1 ]; then
    if [ $UNPAIRED == 0 ]; then
	find $2 -name "${fo}_rRNA*" -print0 | xargs -0 -I {} -P 6 sh -c 'unmerge-paired-reads.sh $0 $1/`basename ${0//.fastq/_1.fq}` $1/`basename ${0//.fastq/_2.fq}`' {} $1
    fi
fi

## keep that as a reminder if that happens again
## sortmerna get confused by the dots as well...
## echo Validating
if [ $UNPAIRED -eq 1 ]; then
    if [ ! -f $1/$fo.fastq ]; then
	mv $1/$fo.* $1/$fo.fq
    fi
fi

## 
echo Gzipping

## compress the output files
find $1 -name "${fo}*.fq" -print0 | xargs -0 -I {} -P 8 gzip -f {}
#printf "%s\0%s" $1/${fo}_1.fq $1/${fo}_2.fq | xargs -0 -I {} -P 2 gzip -f {}

##
echo Done
