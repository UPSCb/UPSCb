#!/usr/bin/env nextflow

params.reads = "$baseDir/*_{1,2}.fq.gz"

def helpMessage() {
	log.info"""
	=========================================
             	 RNAseq pipeline
	=========================================
	
		Example:
	
	nextflow run pl.nf --reads "rnaseq/*_{1,2}.fq.gz" --outDir results --account "Axxxxxx"
	
	Arguments:
		--reads                       Path to input data (must be surrounded with quotes)
		--outDir                      Path to results folder
		--account                     Account for SLURM
		
		
	General options:
		--singleEnd                   Specifies that the input is single end reads
		
		
	FastQC options:
	
	
	sortmeRNA options:
		--smeRNAfastx DEF = true
		--smeRNApairedIn DEF = true
		--smeRNAlog Log details DEF = true

		--smeRNADB Path to databases
	Trimmomatic options:
	trimJar Path to trimmomatic.jar
	trimAdapter Path to adapter file
	trimMode PE (Pair end) or SE (Single end)
	trimPthred  DEF 33
	trimLog = true Record log DEF yes
	trimSeedMismatches Maximum mismatch count DEF = 2
	trimPalClipThreshold  Palindrome Clip Thresgold DEF = 30
	trimSimpleClipThreshold SimpleClipThreshold DEF = 10
	trimWindowSize Number of bases DEF= 5
	trimRequiredQuality Average quality requiered DEF = 20
	trimMinLen minimum length of reads to be kept. DEF = 50
	trimCrop The number of bases to keep, from the start of the read DEF = 0
	trimHeadCrop Number of bases to remove from the start of the read. DEF = 0
	
	Salmon options:
		--salmonIndex                   Path to index
		--salmonMode                    The library type. Default = "IU" 


	""".stripIndent()
}

/**************
***CHECKS

**************/

DUMMY = file('dummy')


if (params.help){
	helpMessage()
	exit 0
}

if (params.aligner != 'star' && params.aligner != 'salmon'){
    exit 1, "Invalid aligner option: ${params.aligner}."
}



/************* END CHECK BLOCK******************/

	Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { read_pairs_ch; read_pairs2_ch	} 
	
/*
if (!params.noSortmeRNA && !params.noTrimmomatic) {
	Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { read_pairs_ch; read_pairs2_ch	} 
	
} else if (params.noSortmeRNA){
	Channel
	.fromFilePairs( params.reads)
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { read_pairs_ch; read_pairs2_ch;sortmerna_unmerged;sortmerna_unmerged2}
} else {
	Channel
	.fromFilePairs( params.reads)
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { read_pairs_ch; read_pairs2_ch;aligner_ch}
}
*/

/*
Step sortmeRNA 
Input: Raw .fastq files
Output: .fastq file without rRNA + FastQC reports + multiqc reports +sortmeRNA
        reports.
*/


if (!params.noSortmeRNA){

	process merger {
		tag "Merge on $pair_id"
		input:
			set pair_id, file(reads) from read_pairs2_ch

		output:
			set pair_id, "${pair_id}_merged.fq" into merged_results
			
		script:	
			//prefix = x.toString() - ~/(_R1)?(_trimmed)?(_merged)?(\.fq)?(\.fastq)?(\.gz)?$/
			"""
			fmerge ${reads} ${pair_id}_merged.fq
			"""
	}

	process sortmerna {
		tag "sortmeRNA on $prefix"
		publishDir "${params.outDir}/${params.step2}", mode: 'copy'
			
		input:
			set pair_id, file(x) from merged_results

		output:
			set pair_id, "*_sortmerna.fq" into sortmerna_results
			file "*_discarded.fq" 
			file "*.log" into sortmerna_logs
			
		script:	
			prefix = x.toString() - ~/(_R1)?(_trimmed)?(_merged)?(\.fq)?(\.fastq)?(\.gz)?$/
			
			def fastx='--fastx'
			def pairedIn = '--paired_in'
			def log = '--log'
			def smrDB = ""
			
			if (!params.smeRNAfastx){	fastx = ''}
			if (!params.smeRNApairedIn){	pairedIn = ''}
			if (!params.smeRNAlog){	log = ''}
			if (params.smeRNADB == "") {smrDB= "\$SORTMERNA_DB"}
		
			"""
			sortmerna --ref $smrDB --reads ${x} \
			--aligned ${prefix}_discarded --other \
			${prefix}_sortmerna $fastx $pairedIn $log -a ${params.cpus}
			"""
	}

	process unmerger {
		tag "Unmerge on $prefix"
		publishDir "${params.outDir}/${params.step2}", mode: 'copy'
		
		input:
			set pair_id, file (x) from sortmerna_results

		output:
			set pair_id, file ("*.fq.gz") into sortmerna_unmerged, sortmerna_unmerged2
			
		script:	
			prefix = x.toString() - ~/(_R1)?(_sortmerna)?(_merged)?(\.fq)?(\.fastq)?(\.gz)?$/
			"""
			funmerge ${x} ${prefix}_sortmerna_1.fq.gz ${prefix}_sortmerna_2.fq.gz
			"""
	}
} else { //end sortmeRNA block
	process no_sortmerna {
		//fakes a sortmerna
		tag "No sortmerna"
		
		input:
			set pair_id, file(reads) from read_pairs2_ch

		output:
			set pair_id, file(reads) into sortmerna_unmerged
		script:
		"""
		"""
	}
}



/**********************
Step3 trimmomatic 
Input: sortmeRNA .fastq files
Output: .fastq file without adapters + FastQC reports 
***********************/
if (!params.noTrimmomatic){
	process trimmomatic {
		publishDir "${params.outDir}/${params.step3}", mode: 'copy'
		input:
			//set file(read1), file(read2) from sortmerna_unmerged
			set pair_id, file(reads) from sortmerna_unmerged

				
		output:
			set file ("*_1.fq.gz"), file ("*_2.fq.gz") into trimmomatic_results
			set pair_id, file ("*_trimmomatic_*.fq.gz") into aligner_ch
			file "*.log" into trimmomatic_logs
			
		script:	
			def read1 =reads[0].toString()
			def read2 =reads[1].toString()
			def prefix = read1.toString() - ~/(_1)?(_sortmerna)?(_2)?(\.fq)?(\.fastq)?(\.gz)?$/
			def mode = params.trimMode
			def trimJar = ""
			def trimAdapter = ""
			def phred = '-phred33'
			def trimLog = "-trimlog ${prefix}_trimmomatic.log"

			def trimSeedMismatches = params.trimSeedMismatches
			def trimPalClipThreshold = params.trimPalClipThreshold
			def trimSimpleClipThreshold = params.trimSimpleClipThreshold
			def trimWindowSize = params.trimWindowSize
			def trimRequiredQuality = params.trimRequiredQuality
			def trimMinLen = params.trimMinLen
			def trimHeadCrop = ''
			def trimCrop = ''
			
			if (params.trimJar == "") {trimJar= "\$TRIMMOMATIC_HOME/trimmomatic.jar"}
			if (params.trimAdapter == "") {trimAdapter= "/build/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa"}
			if (!params.trimLog)
				trimLog=''
			if (params.trimHeadCrop > 0)
				trimHeadCrop= "HEADCROP:${params.trimHeadCrop}"
			if (params.trimCrop > 0)
				trimCrop= "CROP:${params.trimCrop}"

			"""
			//java -jar
			trimmomatic $mode -threads ${params.cpus} $phred $trimLog \
			${read1} ${read2} \
			${prefix}_trimmomatic_1.fq.gz ${prefix}_unpaired_1.fq.gz \
			${prefix}_trimmomatic_2.fq.gz ${prefix}_unpaired_2.fq.gz \
			ILLUMINACLIP:${trimAdapter}:${trimSeedMismatches}:${trimPalClipThreshold}:${trimSimpleClipThreshold} \
			SLIDINGWINDOW:${trimWindowSize}:${trimRequiredQuality} \
			$trimCrop \
			$trimHeadCrop \
			MINLEN:${trimMinLen}
			"""
	}
} else {//end of trimmomatic block
	process no_trimmomatic {
		tag "No trim"
		//fakes a trimmomatic
		
		input:
			set pair_id, file(reads) from sortmerna_unmerged

		output:
			set pair_id, file(reads) into aligner_ch
		script:
		"""
		"""
	}
}



/*
Step4 salmon 
Input: sortmeRNA+trimmomatic .fastq files
Output: count files.
*/

if (params.aligner == 'salmon'){
	process salmon {
		tag "Salmon to ${prefix}"
		publishDir "${params.outDir}/${params.step4}/${pair_id}", mode: 'copy'
		
		input:
			set pair_id, file(reads) from aligner_ch
			//when: !params.noTrimmomatic

		output:
			/*file "quant.sf" into salmon_results*/
			file "*" into aligner_results

		script:	
			def read1 =reads[0].toString()
			def read2 =reads[1].toString()
			//prefix = read1[0].toString() - ~/(_1)?(_sortmerna)?(_trimmomatic)?(_2)?(\.fq)?(\.fastq)?(\.gz)?$/
			"""
			salmon quant -i ${params.salmonIndex} -l ${params.salmonMode} -1 ${read1} \
			-2 ${read2} -p ${params.cpus} --gcBias --output .
			"""
	}  
} else if (params.aligner == 'star' ) {
	process star {
		tag "Star to ${prefix}"
		
		input:
			set pair_id, file(reads) from aligner_ch

		output:
			/*file "quant.sf" into salmon_results*/
			file "*" into aligner_results

		script:	
			def read1 =reads[0].toString()
			def read2 =reads[1].toString()
			//prefix = read1[0].toString() - ~/(_1)?(_sortmerna)?(_trimmomatic)?(_2)?(\.fq)?(\.fastq)?(\.gz)?$/
			"""
			
			//salmon quant -i ${params.salmonIndex} -l ${params.salmonMode} -1 ${read1} \
			//-2 ${read2} -p ${params.cpus} --gcBias --output .
			"""
	}
}

/*****************************************************************************
FASTQC section
Note: inputs and outputs can't be declared dynamically by Nextflow
3 process are needed for each possible FastQC output
*****************************************************************************/

/*
Step FastQC 
Input: Raw .fastq files
Output: FastQC reports + multiqc reports
*/
if (!params.noFastQC){
	process fastqc {

		tag "FASTQC on $pair_id"
		publishDir "${params.outDir}/${params.step1}", mode: 'copy'
		
		input:
			set pair_id, file(reads) from read_pairs_ch

		output:
			file "*_fastqc.{zip,html}" into fastqc_results

		script:	
			def extract='--noextract'
			if (!params.fastqcNoExtract){
			extract = ''
			}
			"""
			fastqc $extract -t ${params.cpus} ${reads}
			"""
	}  

	if (!params.noSortmeRNA){
		process fastqc_sortmerna {
			tag "FASTQC after sortmerna"
			publishDir "${params.outDir}/${params.step2}/${params.fastqcSubDir}", mode: 'copy'
			
			input:
				set file(read1), file(read2) from sortmerna_unmerged2

			output:
				file "*_fastqc.{zip,html}" into fastqc2_results

			script:	
				def extract='--noextract'
				if (!params.fastqcNoExtract){extract = ''}
				"""
				fastqc $extract -t ${params.cpus} ${read1} ${read2}
				"""
		}  
	}
	if (!params.noTrimmomatic){
		process fastqc_trimmomatic {
			tag "FASTQC after trimmomatic"
			publishDir "${params.outDir}/${params.step3}/${params.fastqcSubDir}", mode: 'copy'
			
			input:
				set file(read1), file(read2) from trimmomatic_results

			output:
				file "*_fastqc.{zip,html}" into fastqc3_results

			script:	
				def extract='--noextract'
				if (!params.fastqcNoExtract){extract = ''}
				"""
				fastqc $extract -t ${params.cpus} ${read1} ${read2}
				"""
		}  
	}
} else { //end if(!params.noFastQC)

	//fastqc_results = Channel.create()
	//fastqc2_results = Channel.create()
	//fastqc3_results = Channel.create()
	Channel
	.empty()
	.into { fastqc_results; fastqc2_results; fastqc3_results	} 
}

/*
MultiQC on 
*/
//if (!params.noMultiQC){
	process multiqc{
		tag "$MulqiQc for trimmomatic"
		publishDir "${params.outDir}/MultiQC", mode: 'copy'

		input:
			//file ('*') from multiqc_results.collect()
			file ('*') from aligner_results.collect()
			file ('*') from fastqc_results.collect().ifEmpty(DUMMY)
			//file ('*') from fastqc2_results.collect().ifEmpty(DUMMY)
			//file ('*') from fastqc3_results.ifEmpty(DUMMY)
		
		output:
			file "*"
		
		script:	
			"""
			multiqc . -o ${params.outDir}/MultiQC
			"""
	}
//}


