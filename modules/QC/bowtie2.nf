process FILTERREADS {
	tag "$sampleID"
	label 'bowtie2qc'
	scratch params.scratch

	publishDir "${params.outdir}/qced_fastq", mode: 'copy', pattern: "*.gz", enabled: params.cleanreads
	publishDir "${params.outdir}/bowtie2_logs", mode: 'copy', pattern: "*_log.txt"
	input:
		tuple val(meta), path(reads)
		path(bwt_files)
		each genome
	output:
		tuple val(meta), path("*_clean.fastq.gz"), 		emit: cleaned_reads
		path(bowtie_log), 								emit: filterlog
        path("versions.yml"),          optional: true, 	emit: version

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"
		bowtie_log = sampleID + "_bowtie2_log.txt"
		def input = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]} -U ${reads[2]}"
		"""
		bowtie2 -x $genome \
				$input \
				-S /dev/null \
				--no-unal \
				-p ${task.cpus} \
				--un-gz $unpaired_clean \
				--un-conc-gz ${sampleID}_R%_clean.fastq.gz \
				2> $bowtie_log

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		bowtie2: \$(bowtie2 --version | awk 'FNR==1' |sed 's/.* //')
		END_VERSIONS

		"""
	stub:
		sampleID = meta.id
		bowtie_log = sampleID + "_bowtie2_log.txt"
		if (meta.single_end) {
			"""
			touch ${sampleID}_single_clean.fastq.gz
			touch $bowtie_log
			echo "cleanreads_stub" > versions.yml
			"""
		} else {
							"""
			touch ${sampleID}_R1_clean.fastq.gz
			touch ${sampleID}_R2_clean.fastq.gz
			touch ${sampleID}_single_clean.fastq.gz
			touch $bowtie_log
			echo "filterreads_stub" > versions.yml
			"""
		}
}