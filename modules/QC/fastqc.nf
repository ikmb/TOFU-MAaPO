process FASTQC_raw {

	label 'fastqc'
	tag "$sampleID"
	publishDir "${params.outdir}/FastQC", mode: 'copy', pattern: "*.zip"
	publishDir "${params.outdir}/FastQC", mode: 'copy', pattern: "*.html"

	input:
	tuple val(meta), path(reads)

	output:
	path('*_fastqc.{zip,html}'), emit: fastqc
	//tuple val(meta), path("*_raw.fastq.gz") , emit: reads

	script:
	sampleID = meta.id

	leftnewname = sampleID + "_1_raw.fastq.gz"
    rightnewname = sampleID + "_2_raw.fastq.gz"
	
    if (meta.single_end) {
		"""
        [ ! -f  $leftnewname ] && ln -s ${reads[0]} $leftnewname

		fastqc --quiet --threads ${task.cpus} $leftnewname 

	    """
      } else {
		"""
        [ ! -f  $leftnewname ] && ln -s ${reads[0]} $leftnewname
        [ ! -f  $rightnewname ] && ln -s ${reads[1]} $rightnewname

		fastqc --quiet --threads ${task.cpus} $leftnewname $rightnewname
		"""
	  }
}

process FASTQC_clean {
	tag "$sampleID"
	label 'fastqc'

	publishDir "${params.outdir}/FastQC", mode: 'copy'

	input:
		tuple val(meta), path(reads)

	output:
		path('*_fastqc.{zip,html}'), emit: fastqc

	script:
		sampleID = meta.id

		leftnewname = sampleID + "_1_clean.fastq.gz"
    	rightnewname = sampleID + "_2_clean.fastq.gz"
		unpairednewname = sampleID + "_unpaired_clean.fastq.gz"
		singlenewname = sampleID + "_single_clean.fastq.gz"

		
		if (!meta.single_end) {
		"""
			[ ! -f  $leftnewname ] && ln -s ${reads[0]} $leftnewname
    		[ ! -f  $rightnewname ] && ln -s ${reads[1]} $rightnewname
			[ ! -f  $unpairednewname ] && ln -s ${reads[2]} $unpairednewname

			fastqc --quiet --threads ${task.cpus} $leftnewname $rightnewname $unpairednewname
		"""
		} else {
		"""
			[ ! -f  $singlenewname ] && ln -s ${reads[0]} $singlenewname

			fastqc --quiet --threads ${task.cpus} $singlenewname
		"""
		}
}