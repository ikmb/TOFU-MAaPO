process checkm {

	label 'checkm'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/checkm", mode: 'copy'

	input: 
	    tuple val(sampleID), file(fafile)

	output:
	    file('bins/*')

	shell:
	    """
	    checkm lineage_wf -t ${task.cpus} -x fa . ./bins
	    """
}