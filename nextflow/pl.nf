#!/usr/bin/env nextflow

params.reads = "$baseDir/*_{1,2}.fq.gz"

def helpMessage() {
	log.info"""
	=========================================
             	 RNAseq pipeline
	=========================================
	
		Example:
	
	nextflow run pl.nf --reads "rnaseq/*_{1,2}.fq.gz" --outDir results --account "Axxxxxx"
	
	

	""".stripIndent()
}

params.help = false
if (params.help){
	helpMessage()
	exit 0
}

Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into { read_pairs_ch; read_pairs2_ch }

/*
Step1 FastQC 
Input: Raw .fastq files
Output: FastQC reports + multiqc reports
*/
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
		fastqc $extract -t ${task.cpus} ${reads}
		"""
}  

process multiqc {
	tag "$MultiQC for FastQC after raw data"
	publishDir "${params.outDir}/${params.step1}/${params.multiqcSubDir}", mode: 'copy'

	input:
		file ('*') from fastqc_results.collect()
	
	output:
		file "*"
	
	script:	
		"""
		multiqc .
		"""
}

/*
Step2 sortmeRNA 
Input: Raw .fastq files
Output: .fastq file without rRNA + FastQC reports + multiqc reports +sortmeRNA
        reports.
*/
process merger {
	input:
		set pair_id, file(reads) from read_pairs2_ch

	output:
		file "${pair_id}_merged.fq" into merged_results
		
	script:	
		"""
		fmerge ${reads} ${pair_id}_merged.fq
    """
}

process sortmerna {
	publishDir "${params.outDir}/${params.step2}", mode: 'copy'
		
	input:
		file x from merged_results

	output:
		file "*_sortmerna.fq" into sortmerna_results
		file "*_discarded.fq" 
		file "*.log" into sortmerna_logs
		
	script:	
		prefix = x.toString() - ~/(_R1)?(_trimmed)?(_merged)?(\.fq)?(\.fastq)?(\.gz)?$/
		
		def fastx='--fastx'
		def pairedIn = '--paired_in'
		def log = '--log'
		
		if (!params.smeRNAfastx){		fastx = ''	}
		if (!params.smeRNApairedIn){		pairedIn = ''	}
		if (!params.smeRNAlog){		log = ''	}
  
		"""
		sortmerna --ref ${params.smeRNADB} --reads ${x} \
		--aligned ${prefix}_discarded --other \
		${prefix}_sortmerna $fastx $pairedIn $log -a ${task.cpus}
		"""
}

process multiqc_sortmerna {
	tag "MultiQC for sortmerna logs"
	publishDir "${params.outDir}/${params.step2}/sortmerna_MultiQC", mode: 'copy'

	input:
		file ('*') from sortmerna_logs.collect()
	
	output:
		file "*"
	
	script:	
		"""
		multiqc .
		"""
}

process unmerger {
	publishDir "${params.outDir}/${params.step2}", mode: 'copy'
	
	input:
		file x from sortmerna_results

	output:
		set file ("*_1.fq.gz"), file ("*_2.fq.gz") into sortmerna_unmerged, sortmerna_unmerged2
		
	script:	
		prefix = x.toString() - ~/(_R1)?(_sortmerna)?(_merged)?(\.fq)?(\.fastq)?(\.gz)?$/
		"""
		funmerge ${x} ${prefix}_sortmerna_1.fq.gz ${prefix}_sortmerna_2.fq.gz
		"""
}

process fastqc2 {
	tag "FASTQC after sortmerna"
	publishDir "${params.outDir}/${params.step2}/${params.fastqcSubDir}", mode: 'copy'
	
	input:
		set file(read1), file(read2) from sortmerna_unmerged2

	output:
		file "*_fastqc.{zip,html}" into fastqc2_results

	script:	
		def extract='--noextract'
		if (!params.fastqcNoExtract){
			extract = ''
		}
		"""
		fastqc $extract -t ${task.cpus} ${read1} ${read2}
		"""
}  

process multiqc2 {
	tag "MultiQC for sortmeRNA .fastq files"
	publishDir "${params.outDir}/${params.step2}/${params.multiqcSubDir}", mode: 'copy'

	input:
		file ('*') from fastqc2_results.collect()
	
	output:
		file "*"
	
	script:	
		"""
		multiqc .
		"""
}

/**********************
Step3 trimmomatic 
Input: sortmeRNA .fastq files
Output: .fastq file without adapters + FastQC reports 
***********************/



process trimmomatic {
	publishDir "${params.outDir}/${params.step3}", mode: 'copy'
	input:
		set file(read1), file(read2) from sortmerna_unmerged

	output:
		set file ("*_1.fq.gz"), file ("*_2.fq.gz") into trimmomatic_results
		set file ("*_trimmomatic_1.fq.gz"), file ("*_trimmomatic_2.fq.gz") into trimmomatic_results2
		file "*.log" into trimmomatic_logs
		
	script:	
		def prefix = read1.toString() - ~/(_1)?(_sortmerna)?(_2)?(\.fq)?(\.fastq)?(\.gz)?$/
		def mode = params.trimMode
		def trimAdapter = params.trimAdapter
		def phred = '-phred33'
		def trimLog = "-trimlog ${prefix}_trimmomatic.log"
		def trimJar = params.trimJar
		def trimSeedMismatches = params.trimSeedMismatches
		def trimPalClipThreshold = params.trimPalClipThreshold
		def trimSimpleClipThreshold = params.trimSimpleClipThreshold
		def trimWindowSize = params.trimWindowSize
		def trimRequiredQuality = params.trimRequiredQuality
		def trimMinLen = params.trimMinLen
		def trimHeadCrop = ''
		def trimCrop = ''
		
		if (!params.trimLog)
			trimLog=''
		if (params.trimHeadCrop > 0)
			trimHeadCrop= "HEADCROP:${params.trimHeadCrop}"
		if (params.trimCrop > 0)
			trimCrop= "CROP:${params.trimCrop}"

		"""
		java -jar $trimJar $mode -threads ${task.cpus} $phred $trimLog \
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

/*
DISCARDED TEMPORARY
MultiQC on trimmomatic logs
*/
/*
process multiqc_trimmomatic {
	tag "$MulqiQc for trimmomatic"
	publishDir "${params.outDir}/${params.step3}/trimmomatic_MultiQC", mode: 'copy'

	input:
		file ('*') from trimmomatic_logs.collect()
	
	output:
		file "*"
	
	script:	
		"""
		multiqc .
		"""
}
*/
process fastqc3 {
	tag "FASTQC after trimmomatic"
	publishDir "${params.outDir}/${params.step3}/${params.fastqcSubDir}", mode: 'copy'
	
	input:
		set file(read1), file(read2) from trimmomatic_results

	output:
		file "*_fastqc.{zip,html}" into fastqc3_results

	script:	
		def extract='--noextract'
		if (!params.fastqcNoExtract){
			extract = ''
		}
		"""
		fastqc $extract -t ${task.cpus} ${read1} ${read2}
		"""
}  

process multiqc3 {
	tag "$MultiQC for trimmomatic"
	publishDir "${params.outDir}/${params.step3}/${params.multiqcSubDir}", mode: 'copy'

	input:
		file ('*') from fastqc3_results.collect()
	
	output:
		file "*"
	
	script:	
		"""
		multiqc .
		"""
}

/*
Step4 salmon 
Input: sortmeRNA+trimmomatic .fastq files
Output: count files.
*/
process salmon {
	tag "Salmon to ${prefix}"
	publishDir "${params.outDir}/${params.step4}/${prefix}", mode: 'copy'
	
	input:
		set file(read1), file(read2) from trimmomatic_results2

	output:
		/*file "quant.sf" into salmon_results*/
		file "*"

	script:	
		prefix = read1[0].toString() - ~/(_1)?(_sortmerna)?(_trimmomatic)?(_2)?(\.fq)?(\.fastq)?(\.gz)?$/
		"""
		salmon quant -i ${params.salmonIndex} -l ${params.salmonMode} -1 ${read1} -2 ${read2} -p ${task.cpus} --output .
		"""
}  


