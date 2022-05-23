process FILTERREADS_SE {

	label 'bowtie2'
	scratch params.scratch

	if(params.cleanreads){
		publishDir "${params.outdir}/qced_files", mode: 'copy'
	}

	input:
		tuple val(meta), path(reads)
		path(bwt_files)
		each genome
	output:
		tuple val(meta), path(unpaired_clean), emit: cleaned_reads
		path(bowtie_log), emit: filterlog

	script:
		sampleID = meta.id

		//left_clean = sampleID + "_R1_clean.fastq.gz"
		//right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"
		bowtie_log = sampleID + "_bowtie2_log.txt"
		
		"""
		bowtie2 --met-stderr \
				-x $genome \
				-U $reads \
				-S /dev/null \
				--no-unal \
				-p ${task.cpus} \
				--un-gz $unpaired_clean \
				--un-conc-gz ${sampleID}_R%_clean.fastq.gz \
				2> $bowtie_log
		"""
}

process FILTERREADS_PE {

	label 'bowtie2'
	scratch params.scratch

	if(params.cleanreads){
		publishDir "${params.outdir}/qced_files", mode: 'copy'
	}

	input:
		tuple val(meta), path(reads), path(unpaired)
		path(bwt_files)
		each genome
	output:
		tuple val(meta), path("*_clean.fastq.gz"), emit: cleaned_reads
		path(bowtie_log), emit: filterlog

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"
		bowtie_log = sampleID + "_bowtie2_log.txt"

		"""
		bowtie2 --met-stderr \
				-x $genome \
				-1 ${reads[0]} \
				-2 ${reads[1]} \
				-U ${unpaired} \
				-S /dev/null \
				--no-unal \
				-p ${task.cpus} \
				--un-gz $unpaired_clean \
				--un-conc-gz ${sampleID}_R%_clean.fastq.gz \
				2> $bowtie_log
		"""
}