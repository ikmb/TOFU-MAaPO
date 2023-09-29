
process FASTP {
	tag "$sampleID"
	label 'fastp'
	scratch params.scratch
	publishDir "${params.outdir}/FASTP", mode: 'copy', pattern: "*.html"

	input:
		tuple val(meta), path(reads)
		
	output:
		
		tuple val(meta), path("*_trimmed.fastq.gz"), emit: filterReads
		path("versions.yml"), emit: version
		path("*fastp.{html,json}"), emit: html_report
	script:
		sampleID = meta.id

		html_report = sampleID + ".fastp.html"
		json_report = sampleID + ".fastp.json"
		log_report = sampleID + ".fastp.log"
		leftnewname = sampleID + "_1_raw.fastq.gz"
		rightnewname = sampleID + "_2_raw.fastq.gz"

		left_trimmed = sampleID + "_1_trimmed.fastq.gz"
		right_trimmed = sampleID + "_2_trimmed.fastq.gz"
		unpaired = sampleID + "_unpaired_trimmed.fastq.gz"

		if (meta.single_end) {
			"""
			[ ! -f  $leftnewname ] && ln -s ${reads} $leftnewname

			fastp \
				--thread=${task.cpus} \
				--in1=${leftnewname} \
				--out1=${left_trimmed} \
				--json=${json_report} \
				--html=${html_report} 2> ${log_report}

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			fastp: \$(fastp --version 2>&1 | awk '{print \$2}' )
			END_VERSIONS
			"""
		} else {
			"""
			[ ! -f  $leftnewname ] && ln -s ${reads[0]} $leftnewname
			[ ! -f  $rightnewname ] && ln -s ${reads[1]} $rightnewname

			fastp \
				--thread=${task.cpus} \
				--in1=${leftnewname} \
				--in2=${rightnewname} \
				--out1=${left_trimmed} \
				--out2=${right_trimmed} \
				--html=${html_report} \
				--json=${json_report} \
				--unpaired1 ${sampleID}_unpaired1.fastq.gz \
				--unpaired2 ${sampleID}_unpaired2.fastq.gz \
				--detect_adapter_for_pe

			cat ${sampleID}_unpaired1.fastq.gz ${sampleID}_unpaired2.fastq.gz > ${unpaired}

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			fastp: \$(fastp --version 2>&1 | awk '{print \$2}' )
			END_VERSIONS
			"""
		}
}