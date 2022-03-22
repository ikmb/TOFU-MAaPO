/*
# These are processes to keep the cardinality and name structure the same
# with the workflow that removes host reads
*/

process COLLECTOR_PE {

	scratch params.scratch

	if(params.cleanreads){
		publishDir "${params.outdir}/reads_clean", mode: 'copy'
	}

	input:
		tuple val(meta), path(reads), path(unpaired)
    output:
		tuple val(meta), path("*_clean.fastq.gz"), emit: cleaned_reads

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		"""
        ln -s ${reads[0]} $left_clean
    	ln -s ${reads[1]} $right_clean
    	ln -s ${unpaired} $single_clean
	    """
}

process COLLECTOR_SE {

	scratch params.scratch

	if(params.cleanreads){
		publishDir "${params.outdir}/reads_clean", mode: 'copy'
	}

	input:
		tuple val(meta), path(reads)
    output:
		tuple val(meta), path("*_clean.fastq.gz"), emit: cleaned_reads

	script:
		sampleID = meta.id

		//left_clean = sampleID + "_R1_clean.fastq.gz"
		//right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

        """
        ln -s ${reads[0]} $single_clean
        """
}