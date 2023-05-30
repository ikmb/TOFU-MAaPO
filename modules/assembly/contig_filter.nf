process FILTERCONTIGS {
	scratch params.scratch
	tag "$sampleID"
	label 'default'

	input:
    	tuple val(meta), file(finalcontigs), path(reads)

	output:
	    tuple val(meta), file(fcontigs_filtered), path(reads), 	emit: contigs, 			optional: true
	    path(fcontigs_filtered), 								emit: filteredcontig,	optional: true
	    tuple val(meta), file(fcontigs_filtered), 				emit: magscot_contigs, 	optional: true
        path("versions.yml"),  									emit: versions,			optional: true

	script:
		sampleID = meta.id
		fcontigs_filtered = sampleID + '_fcontigsfiltered.fa'

		"""
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/contigfilterbylen.py ${params.contigsminlength} $finalcontigs > intermediate_file.fa

		if [ -s intermediate_file.fa ]; then
        	# The file is not-empty.
        	mv intermediate_file.fa $fcontigs_filtered
		else
        	# The file is empty.
			echo "file is empty"
		fi

		cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        Python: \$(/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 --version | sed -e "s/Python //g" )
        END_VERSIONS
		"""
}