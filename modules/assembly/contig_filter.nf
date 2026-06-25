process FILTERCONTIGS {
	scratch params.scratch
	tag "$coassemblygroup"
	label 'seqkit'
	label 'short_run'
	input:
		tuple val(coassemblygroup), file(finalcontigs)

	output:
		tuple val(coassemblygroup), file(fcontigs_filtered),	emit: contigs, 			optional: true
		path(fcontigs_filtered), 								emit: filteredcontig,	optional: true
		tuple val(coassemblygroup), file(fcontigs_filtered), 	emit: magscot_contigs, 	optional: true
		path("versions.yml"),  									emit: versions,			optional: true

	script:

		fcontigs_filtered = coassemblygroup + '_fcontigsfiltered.fa'

		"""
		seqkit seq -m ${params.contigsminlength} $finalcontigs > intermediate_file.fa

		if [ -s intermediate_file.fa ]; then
			# The file is not-empty.
			mv intermediate_file.fa $fcontigs_filtered
		else
			# The file is empty.
			echo "file is empty"
		fi

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		SeqKit: \$(seqkit version | sed -e "s/^seqkit //")
		END_VERSIONS
		
		"""
	stub:
		fcontigs_filtered = coassemblygroup + '_fcontigsfiltered.fa'

		"""
		touch $fcontigs_filtered
		echo "FILTERCONTIGS_stub" > versions.yml
		"""
}
