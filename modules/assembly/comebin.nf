process COMEBIN {

	label 'comebin'
	label 'gpu'
	scratch params.scratch
	tag { meta.id }
	cache 'lenient'
	errorStrategy { task.exitStatus == 42 ? 'ignore' : task.attempt <= 2 ? 'retry' : 'finish' } //this would not publish failed logs

	publishDir "${params.outdir}/comebin", mode: 'copy', enabled: params.publish_rawbins,
        saveAs: { filename -> "${meta.id}/${filename}" }

	publishDir "${params.outdir}/comebin",
        mode: 'copy',
        pattern: "*_comebin_failed.txt",
        saveAs: { filename -> "${meta.id}/${filename}" }

	container  { params.gpu ?	"docker://eikematthias/comebin:1.0.4" : 
								"docker://eikematthias/comebin:1.0.4" } //TODO: create small non-gpu container
	containerOptions { params.gpu ? '--nv' : '' }

	input: 
		tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)

	output:
		tuple val(meta), file("${meta.id}_bin.*.fa"), optional: true, emit: comebinout
		tuple val(meta), file("${meta.id}_comebin_magscot_contigs_to_bin.tsv"), optional: true, emit: magscot_contigbinlist
		tuple val(meta), file("${sampleID}_comebin_failed.txt"), optional: true, emit: comebine_failed
		path("versions.yml"),          optional: true, emit: versions

	script:
		def sampleID = meta.id
		def comebin_contigs_to_bin = sampleID + '_comebin_output/comebin_res/comebin_res.tsv'
		def formatted_contigs_to_bin = sampleID + '_comebin_magscot_contigs_to_bin.tsv'
		//def engine = params.gpu ? 'gpu' : 'cpu'
			"""
			if run_comebin.sh \
				-a $fcontigs \
				-p . \
				-o ${sampleID}_comebin_output \
				-t ${task.cpus} \
				-n 6 > >(tee ${sampleID}_comebin.log) 2>&1
			then
    		:
			else
			    comebin_exit=\$?

				if grep -Fq "RuntimeError: selected index k out of range" ${sampleID}_comebin.log; then
					{
						echo "WARNING: COMEBin skipped sample ${sampleID}: insufficient contigs for network training."
						echo
						cat ${sampleID}_comebin.log
					} > ${sampleID}_comebin_failed.txt

					exit 0
				fi

				exit "\$comebin_exit"
			fi

			# old implementation
			if [[ ! -s "$comebin_contigs_to_bin" ]]; then
				echo "WARNING: COMEBin produced no usable bins for ${sampleID}; continuing without COMEBin results." >&2
				exit 42
			fi

			find ${sampleID}_comebin_output -type f -name '*.fa' -print0 | while IFS= read -r -d '' bin; do
				dest="${sampleID}_bin.\$(basename "\${bin}")"
				ln -sfn "\${bin}" "\${dest}"
			done
			awk 'NR>1 {print "${sampleID}_comebin_"\$2"\t"\$1"\tcomebin"}' $comebin_contigs_to_bin > $formatted_contigs_to_bin

			cat <<-END_VERSIONS> versions.yml
			"${task.process}":
			comebin: 1.0.4
			END_VERSIONS
			"""
	stub:
		def sampleID = meta.id
		def comebin_contigs_to_bin = sampleID + '_comebin_contigs_to_bin.tsv'
		def formatted_contigs_to_bin = sampleID + '_comebin_magscot_contigs_to_bin.tsv'

		"""
		touch $comebin_contigs_to_bin
		touch ${sampleID}_bin.stub.fa
		touch $formatted_contigs_to_bin

		echo "comebin_stub" > versions.yml
	"""
}
