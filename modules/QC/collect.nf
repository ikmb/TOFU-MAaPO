/*
 * This is a process to keep the cardinality and name structure the same
 * to the workflow that removes host reads
 */

process COLLECTOR {

	label 'default'
	scratch params.scratch
	tag "$sampleID"

	if(params.cleanreads){
		publishDir "${params.outdir}/reads_clean", mode: 'copy'
	}

	input:
		tuple val(meta), path(reads)
    output:
		tuple val(meta), path("*_clean.fastq.gz"), emit: cleaned_reads

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		"""
        ln -s ${reads[0]} $left_clean
    	if [ -f "${reads[1]}" ]; then ln -s ${reads[1]} $right_clean; fi
    	if [ -f "${reads[2]}" ]; then ln -s ${reads[2]} $single_clean; fi
	    """
}
