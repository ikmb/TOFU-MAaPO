process checkm {

	label 'checkm'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/checkm/${sampleID}", mode: 'copy'

	input: 
	    tuple val(meta), file(fafile)

	output:
	    file('bins/*')
		file(outputtable)
	shell:
		sampleID = meta.id
		outputtable = sampleID + "_checkm_table.tsv"
	    """
	    checkm lineage_wf -t ${task.cpus} -x fa . ./bins --tab_table --file ${outputtable}
	    """
}