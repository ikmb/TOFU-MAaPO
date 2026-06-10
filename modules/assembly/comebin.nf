process COMEBIN {

	label 'comebin'
	label 'long_run'
	scratch params.scratch
	tag "$sampleID"
	cache 'lenient'

	container  "docker://eikematthias/comebin:1.0.4"

	publishDir "${params.outdir}/comebin", mode: 'copy', enabled: params.publish_rawbins,
        saveAs: { filename -> "${meta.id}/${filename}" }

	input: 
		tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)

	output:
		tuple val(meta), file("${sampleID}_bin.*.fa"), optional: true, emit: comebinout
		tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id
		comebin_contigs_to_bin = sampleID + '_comebin_output/comebin_res/comebin_res.tsv'
		formatted_contigs_to_bin = sampleID + '_comebin_magscot_contigs_to_bin.tsv'
		//def engine = params.gpu ? 'gpu' : 'cpu'
			"""
			run_comebin.sh \
				-a $fcontigs \
				-p . \
				-o ${sampleID}_comebin_output \
				-t ${task.cpus} \
				-n 6 

			awk 'NR>1 {print "${sampleID}_comebin_"\$2"\t"\$1"\tcomebin"}' $comebin_contigs_to_bin > $formatted_contigs_to_bin

			cat <<-END_VERSIONS> versions.yml
			"${task.process}":
			comebin: 1.0.4
			END_VERSIONS
			"""
	stub:
		sampleID = meta.id
		bed_file = sampleID + '.bed'
		comebin_contigs_to_bin = sampleID + '_comebin_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_comebin_magscot_contigs_to_bin.tsv'

		"""
		touch $comebin_contigs_to_bin
		touch $formatted_contigs_to_bin

		echo "comebin_stub" > versions.yml
		"""
}
