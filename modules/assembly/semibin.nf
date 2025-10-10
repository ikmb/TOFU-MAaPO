process SEMIBIN {

	label 'semibin'
	label 'very_long_run'
	scratch params.scratch
	tag "$sampleID"
	cache 'lenient'
	publishDir "${params.outdir}/semibin/${sampleID}", mode: 'copy', enabled: params.publish_rawbins

	input: 
		tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)

	output:
		tuple val(meta), file("${sampleID}_bin.*.fa"), optional: true, emit: semibinout
		tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id
		semibin_contigs_to_bin = sampleID + '_semibin_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_semibin_magscot_contigs_to_bin.tsv'
			"""
			SemiBin2 single_easy_bin \
				-i $fcontigs \
				-b $mappingbam \
				-o ${sampleID}_semibin_output \
				--environment ${params.semibin_environment} 

			set +e
			{
				grep '>' ${sampleID}_semibin_output/output_recluster_bins/*.fa | cut -d '/' -f 3 | sed 's/:>/\\ /' > $semibin_contigs_to_bin
			}
			set -e

			awk '{print \$1"\t"\$2"\tsemibin2"}' $semibin_contigs_to_bin > $formatted_contigs_to_bin

			cat <<-END_VERSIONS> versions.yml
			"${task.process}":
			SemiBin: \$(SemiBin --version 2>&1)
			END_VERSIONS

			"""

}