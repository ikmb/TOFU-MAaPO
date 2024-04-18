
process DREP {

	label 'drep'
	scratch params.scratch
	//tag "${sampleID}_${meta.assembler}"
	publishDir "${params.outdir}/drep", mode: 'copy'

	input:
		tuple val(coassemblygroup), val(meta), file(formatted_contigs_to_bin), file(samplehmm), file(fcontigs_filtered)
	output:
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered), emit: drep_results
		path("versions.yml"),          optional: true, emit: versions
	script:
		//sampleID = meta.id
		refined_contigs_to_bins = sampleID + '.refined.contig_to_bin.out'
		stats_outfile = sampleID + '.refined.out'
		full_stats = sampleID + '.scores.out'
	"""
		#drep

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python3 --version | sed -e "s/Python //g" )
		END_VERSIONS
	"""
}
