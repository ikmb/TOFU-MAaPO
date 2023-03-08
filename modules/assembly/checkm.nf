process checkm {

	label 'checkm'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/checkm/${sampleID}", mode: 'copy'

	input: 
	    tuple val(meta), file(fafile)

	output:
	    file('bins/*')

	shell:
		sampleID = meta.id
	    """
	    checkm lineage_wf -t ${task.cpus} -x fa . ./bins
	    """
}