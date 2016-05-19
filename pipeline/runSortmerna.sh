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
KEEP=1
useMtSSU=0
UNPAIRED=0
PROC=16
DBS=

## usage
usage(){
echo >&2 \
"
	Usage: runSortmerna.sh [option] <out dir> <tmp dir> <forward fastq.gz> <reverse fastq.gz>
	
	Options:
                -d define your dbs (semi-colon separated)
                -k drop the rRNA (only for v1.9, default to keep them)
                -m run against mtSSU in addition (only for v1.9)
                -p number of threads to be used (default $PROC)
                -u single end data (in that case only the forward fastq is needed)

         Note:
               1) The SORTMERNADIR environment variable needs to be set
               2) Only SortMeRna version 1.9 and 2.x are supported (2.x is default)
               3) -m is not applicable if -d is set
"
	exit 1
}

## load the module
module load bioinfo-tools

## Does not work on uppmax - umea has an empty result
## while uppmax is verbose.
## avail=$( module avail sortmerna 2>&1 > /dev/null)
## avail=`echo $avail | tr -d [:blank:]`
## if [ ! -z $avail ]; then
##  module load sortmerna
##  sortmerna --version
##fi

## record the SORTMERNADIR if it exists
STOREENV=
if [ ! -z $SORTMERNADIR ]; then
  STOREENV=$SORTMERNADIR
fi

## try to load or echo
module load sortmerna || {
  echo "No sortmerna as module"

  ## then check for availability
  tool=`which sortmerna 2>/dev/null`
  if [ ! -z $tool ] && [ -f $tool ] && [ -x $tool ]; then
    echo "sortmerna available"
  else
    echo "ERROR: INSTALL SortMeRna"
    usage
  fi
}

# restore the env if it existed
if [ ! -z $STOREENV ]; then
  export SORTMERNADIR=$STOREENV
fi

## check for sortmerna version
is1dot9=`sortmerna --version 2>&1 | grep version | grep 1.9 | wc -c`
is2dotx=`sortmerna --version 2>&1 | grep "version 2." | wc -c`

if [ $is1dot9 == 0 ] && [ $is2dotx  == 0 ]; then
  echo "Only version 1.9 and 2.x are supported"
  usage
fi

## get the options
while getopts d:kmp:u option
do
        case "$option" in
      d) DBS=$OPTARG;;
	    k) KEEP=0;;
	    m) useMtSSU=1;;
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

## set the default dbs
if [ ! -z $DBS ]; then
  dbs=${DBS//;/ }
  dbNum=`echo $DBS | awk -F";" '{print NF}'`
else
  if [ $is2dotx != 0 ]; then
    db5s=$SORTMERNADIR/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNADIR/automata/rfam-5s-database-id98
    db58s=$SORTMERNADIR/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNADIR/automata/rfam-5.8s-database-id98
    db16sa=$SORTMERNADIR/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNADIR/automata/silva-arc-16s-database-id95
    db16s=$SORTMERNADIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNADIR/automata/silva-bac-16s-database-id90
    db18s=$SORTMERNADIR/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNADIR/automata/silva-euk-18s-database-id95
    db23sa=$SORTMERNADIR/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNADIR/automata/silva-arc-23s-database-id98
    db23s=$SORTMERNADIR/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNADIR/automata/silva-bac-23s-database-id98
    db28s=$SORTMERNADIR/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNADIR/automata/silva-euk-28s-database-id98
    dbs="$db5s:$db58s:$db16sa:$db16s:$db18s:$db23sa:$db23s:$db28s"
  #if [ ! -f $SORTMERNADIR/automata/rfam-5s-database-id98.stats ]; then
  #  echo "No indexes found, creating indexes in folder $SORTMERNADIR/automata"
  #  indexdb_rna --ref $dbs
  #fi
  else
    db5s=$SORTMERNADIR/rRNA_databases/rfam-5s-database-id98.fasta
    db58s=$SORTMERNADIR/rRNA_databases/rfam-5.8s-database-id98.fasta
    db16sa=$SORTMERNADIR/rRNA_databases/silva-arc-16s-database-id95.fasta
    db16s=$SORTMERNADIR/rRNA_databases/silva-bac-16s-database-id85.fasta
    db18s=$SORTMERNADIR/rRNA_databases/silva-euk-18s-database-id95.fasta
    db23sa=$SORTMERNADIR/rRNA_databases/silva-arc-23s-database-id98.fasta
    db23s=$SORTMERNADIR/rRNA_databases/silva-bac-23s-database-id98.fasta
    db28s=$SORTMERNADIR/rRNA_databases/silva-euk-28s-database-id98.fasta
    dbNum=8
    dbs="$db5s $db58s $db16sa $db16s $db18s $db23sa $db23s $db28s"
  fi

  ## Add the mtSSU
  if [ $is1dot9 != 0 ] && [ $useMtSSU == 1 ]; then
    mtSSU=$SORTMERNADIR/rRNA_databases/mtSSU_UCLUST-95-identity.fasta
    dbs="$dbs $mtSSU"
    dbNum=9
  fi
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
  merge-paired-reads.sh $2/$f1 $2/$f2 $2/$fm
fi

##
if [ $UNPAIRED == 0 ]; then
    echo Pre-cleaning
    rm -f $2/$f1 $2/$f2
else
    echo "TODO: Cleaning needs implementing for single end sequencing"
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
opt="-a $PROC"

if [ $KEEP == 1 ] && [ $is1dot9 != 0 ]; then
  opt="$opt --bydbs --accept $2/${fo}_rRNA"
fi 

## run
if [ $UNPAIRED == 0 ]; then
  if [ $is2dotx != 0 ]; then
    sortmerna --ref $dbs --reads $2/$fm --other $2/$fo --log --paired_in --fastx $opt --sam --num_alignments 1 --aligned $2/${fo}_rRNA
  else
    sortmerna -n $dbNum --db $dbs --I $2/$fm --other $2/$fo --log $1/$fo --paired-in $opt
  fi  
else
  if [ $is2dotx != 0 ]; then
    sortmerna --ref $dbs --reads $2/$f1 --other $1/$fo --log $opt --sam --fastx --num_alignments 1 --aligned $2/${fo}_rRNA
  else
    sortmerna -n $dbNum --db $dbs --I $2/$f1 --other $1/$fo --log $1/$fo $opt
  fi
fi

## deinterleave it
if [ $UNPAIRED == 0 ]; then
    ## sortmerna get confused by dots in the filenames
    if [ ! -f $2/$fo.fastq ]; then
	    mv $2/$fo.* $2/$fo.fastq
    fi
    unmerge-paired-reads.sh $2/$fo.fastq $1/${fo}_1.fq $1/${fo}_2.fq
fi

## cleanup
echo Post-Cleaning

if [ $is2dotx != 0 ]; then
  ## mv the rRNA, fastq and log back
  mv $2/${fo}_rRNA.* $1
fi

## rm the tmp
if [ $UNPAIRED == 0 ]; then
    rm -f $2/$fm $2/$fo.fastq
else
    rm -f $2/$f1
fi

## deinterleave the rest if needed
if [ $KEEP == 1 ]; then
    if [ $UNPAIRED == 0 ]; then
	find $2 -name "${fo}_rRNA*" -print0 | xargs -0 -I {} -P 6 sh -c 'unmerge-paired-reads.sh $0 $1/`basename ${0//.fastq/_1.fq}` $1/`basename ${0//.fastq/_2.fq}`' {} $1
    fi
fi

## keep that as a reminder if that happens again
## sortmerna get confused by the dots as well...
## echo Validating
if [ $UNPAIRED == 1 ]; then
    if [ ! -f $1/$fo.fastq ]; then
      mv $1/$fo.* $1/$fo.fq
    fi
fi

## 
echo Gzipping

## compress the output files
find $1 -name "${fo}*.fq" -print0 | xargs -0 -I {} -P 8 gzip -f {}
if [ $is2dotx != 0 ]; then
  find $1 -name "${fo}_rRNA.[f,s]*" -print0 | xargs -0 -I {} -P 8 gzip -f {}
fi

## TODO if unpaired, then the input file is still called fastq and hence not compressed
## FIXME

#printf "%s\0%s" $1/${fo}_1.fq $1/${fo}_2.fq | xargs -0 -I {} -P 2 gzip -f {}

##
echo Done
