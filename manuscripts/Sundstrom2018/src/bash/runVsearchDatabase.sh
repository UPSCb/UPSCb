#!/bin/bash -l

#SBATCH -p node
#SBATCH -n 20
#SBATCH --mem 64000
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=ALL

##be verbose and stop on error
set -ex

##usage
usage(){
    echo >&2 \
"Usage: $0 [options] <output folder> <file> [file] ...

Submit all samples you want to concatenate to construct a database of the OTUs throughout all samples.

Options:
       -t  type of input data (16S or ITS, default: 16S)
       -R  reference based chimera checking ACTIVATED (if no special database desired, skip -r, defaults are set)
       -r  specify database file for reference based chimera checking (defaults for 16S and ITS are\ set)
       -i  sequence identity (if not specified, it will be set to the default value of 0.97 for 16S and for ITS)
       -d  specify database for taxonomy assignment (defaults are Silva 123 for 16S and Unite V9 for ITS)
       -m  specify different mapping file than default linked to Silva/unite, has to be specified if -d is

Note:
You need to set the UPSCb env. variable to your UPSCb git checkout directory
"
    exit 1
}

##load modules
module load bioinfo-tools vsearch Qiime blast/2.2.26

##defaults
typ=16S
def=/mnt/picea/projects/metaseq/asjoedin/qiime_databases/usearch_chimera_bacteria/rdp_gold.fa
ref=
ident=0.97
daba=/mnt/picea/storage/reference/Qiime/silva/Silva_128/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta
mapfile=/mnt/picea/storage/reference/Qiime/silva/Silva_128/SILVA_128_QIIME_release/taxonomy/16S_only/97/majority_taxonomy_7_levels.txt

##options
while getopts t:Rr:i:d:m: opt;
do
    case $opt in
	t) typ=$OPTARG;;
	R) ref=$def;;
	r) ref=$OPTARG;;
	i) ident=$OPTARG;;
	d) daba=$OPTARG;;
	m) mapfile=$OPTARG;;
	\?) usage;;
    esac
done

shift `expr $OPTIND - 1`

#arguments
if [ $# -le 2 ]; then
    echo "This function takes at least two arguments, the output folder and at least one valid fasta file"
    usage
fi

if [ ! -f $2 ]; then
    echo "The last argument(s) must be (a) valid fasta file(s)"
    usage
fi



##create variables depending on input options
if [ $typ == "ITS" ]; then
    if [ $daba == "/mnt/picea/storage/reference/Qiime/silva/Silva_128/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta" ]; then
	daba=/mnt/picea/storage/reference/Qiime/unite/UNITE_7.2/sh_refs_qiime_ver7_dynamic_s_28.06.2017.fasta
	mapfile=/mnt/picea/storage/reference/Qiime/unite/UNITE_7.2/sh_taxonomy_qiime_ver7_dynamic_s_28.06.2017.txt
    fi    
    if [ $ref == "/mnt/picea/projects/metaseq/asjoedin/qiime_databases/usearch_chimera_bacteria/rdp_gold.fa" ]; then
	ref=/mnt/picea/projects/metaseq/asjoedin/qiime_databases/usearch_chimera_fungi/uchime_reference_dataset_01.01.2016/ITS1_ITS2_datasets/uchime_sh_refs_dynamic_develop_985_01.01.2016.ITS1.fasta
    fi
fi

DIR=$1
shift

##concatenate all input fasta files

cat $@ >> $DIR/Pooled_Samples.fa

PS=$DIR/Pooled_Samples.fa

##sort by length
vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --sortbylength $PS --output $DIR/sortsiz.fa

rm -f $PS

##dereplicate
vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --derep_fulllength $DIR/sortsiz.fa --output $DIR/derep.fa --sizeout

rm -f $DIR/sortsiz.fa

##sort by cluster size
vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --sortbysize $DIR/derep.fa --output $DIR/sorted.fa --minsize 2 --sizein --sizeout

rm -f $DIR/derep.fa

##cluster by identity
vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --cluster_size $DIR/sorted.fa --id $ident --strand both --centroids $DIR/centroids.fa --relabel_sha1 --sizein --sizeout

rm -f $DIR/sorted.fa

##filter out chimeras globally
vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --uchime_denovo $DIR/centroids.fa --nonchimeras $DIR/nonchimeras.fa --chimeras $DIR/chimeras.fa --sizein --sizeout

#database chimera checking (optional)
if [ $ref ]; then
    vsearch --threads=$SLURM_JOB_CPUS_PER_NODE --uchime_ref $DIR/nonchimeras.fa --db $ref --nonchimeras $DIR/refnonchimeras.fa --chimeras $1/refchimeras.fa --sizein --sizeout
fi

##Assign taxonomy to all centroid sequences
parallel_assign_taxonomy_blast.py -i $DIR/nonchimeras.fa -o $DIR/taxonomy -T --jobs_to_start 10 --reference_seqs_fp $daba --id_to_taxonomy_fp $mapfile
