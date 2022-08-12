process filtercontigs {
	scratch params.scratch
	tag "$sampleID"

	input:
    	tuple val(sampleID), file(finalcontigs), path(reads)

	output:
	    tuple val(sampleID), file(fcontigs_filtered), path(reads), emit: contigs
	    path(fcontigs_filtered), emit: filteredcontig

	    tuple val(sampleID), file(fcontigs_filtered), emit: magscot_contigs

	script:
		fcontigs_filtered = sampleID + '_fcontigsfiltered.fa'

		"""
		python3 ${baseDir}/bin/contigfilterbylen.py ${params.contigsminlength} $finalcontigs > $fcontigs_filtered
		"""
}