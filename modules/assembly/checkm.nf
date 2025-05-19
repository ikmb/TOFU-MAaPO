process checkm {

	label 'checkm'
	label 'long_run'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/checkm/${sampleID}", mode: 'copy'

	input: 
		tuple val(meta), file(fafile)

	output:
		path('bins/*'),         optional: true, emit: alloutputs
		path(outputtable), 		optional: true, emit: checkmtable
		path("versions.yml"),	optional: true, emit: versions
	shell:
		sampleID = meta.id
		outputtable = sampleID + "_checkm_table.tsv"
		"""
		checkm lineage_wf -t ${task.cpus} -x fa . ./bins --tab_table --file ${outputtable}

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		checkm: \$(checkm -h | awk 'NR==2{print \$3}' | sed -e 's/v//g')
		END_VERSIONS
		"""
}