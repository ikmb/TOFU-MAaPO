process filtercontigs {
	scratch params.scratch
	tag "$sampleID"
	label 'default'

	input:
    	tuple val(meta), file(finalcontigs), path(reads)

	output:
	    tuple val(meta), file(fcontigs_filtered), path(reads), emit: contigs
	    path(fcontigs_filtered), emit: filteredcontig

	    tuple val(meta), file(fcontigs_filtered), emit: magscot_contigs

	script:
		sampleID = meta.id

		fcontigs_filtered = sampleID + '_fcontigsfiltered.fa'

		"""
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/contigfilterbylen.py ${params.contigsminlength} $finalcontigs > $fcontigs_filtered
		"""
}