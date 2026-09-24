process SEMIBIN {

	label 'semibin'
	label 'very_long_run'
	scratch params.scratch
	tag "${meta.id}"
	cache 'lenient'
	publishDir "${params.outdir}/semibin", mode: 'copy', enabled: params.publish_rawbins,
        saveAs: { filename -> "${meta.id}/${filename}" }

	publishDir "${params.outdir}/semibin",
		mode: 'copy',
		pattern: "*_semibin_failed.txt",
		saveAs: { filename -> "${meta.id}/${filename}" }
	label 'gpu'
	container  { params.gpu ?	"docker://eikematthias/semibin_gpu:2.2.0" : 
								"docker://quay.io/biocontainers/semibin:2.2.0--pyhdfd78af_0" }
	containerOptions { params.gpu ? '--nv' : '' }


	input: 
		tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)

	output:
		tuple val(meta), file("${sampleID}_bin.*.fa"), optional: true, emit: semibinout
		tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
		tuple val(meta), file("${sampleID}_semibin_failed.txt"), optional: true, emit: semibin_failed
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id
		semibin_contigs_to_bin = sampleID + '_semibin_output/contig_bins.tsv'
		formatted_contigs_to_bin = sampleID + '_semibin_magscot_contigs_to_bin.tsv'
		def engine = params.gpu ? 'gpu' : 'cpu'
		"""
		SEMIBIN_LOG="${sampleID}_semibin.log"
		SEMIBIN_FAILED="${sampleID}_semibin_failed.txt"

		if SemiBin2 single_easy_bin \
			-i $fcontigs \
			-b $mappingbam \
			-o ${sampleID}_semibin_output \
			--engine ${engine} \
			-t ${task.cpus} \
			--environment ${params.semibin_environment} 
		then
			:
		else
			semibin_status=\$?

			# SemiBin execution log
			cat ${sampleID}_semibin_output/SemiBinRun.log > "\$SEMIBIN_LOG"
			cat "\$SEMIBIN_LOG" >&2
			if grep -Fq \
			    -e "all are shorter than 2500 basepairs" \
    			-e "contain(s) at least 1000 basepairs" \
				"\$SEMIBIN_LOG"; then

				echo "WARNING: SemiBin skipped sample ${sampleID}: all contigs are shorter than 2500 bp." > "\$SEMIBIN_FAILED"

				# Also report it in the Nextflow log
				cat "\$SEMIBIN_FAILED" >&2

				# A single sample should not kill the pipeline
				exit 0
			fi
			# Any OTHER SemiBin error is still a real error
			exit "\$semibin_status"
		fi

		awk 'NR>1 {print "semibin2_"\$2"\t"\$1"\tsemibin2"}' $semibin_contigs_to_bin > $formatted_contigs_to_bin

		cat <<-END_VERSIONS> versions.yml
		"${task.process}":
		SemiBin2: \$(SemiBin2 --version 2>&1)
		END_VERSIONS
		"""
	stub:
		sampleID = meta.id
		semibin_contigs_to_bin = sampleID + '_semibin_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_semibin_magscot_contigs_to_bin.tsv'
		"""
		touch $semibin_contigs_to_bin
		touch $formatted_contigs_to_bin
		echo "SEMIBIN_stub" > versions.yml
		"""

}
