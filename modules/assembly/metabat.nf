process METABAT {

    label 'metabat2'
    scratch params.scratch
    tag "$sampleID"
    publishDir "${params.outdir}/${sampleID}/Metabat2", mode: 'copy'

    input: 
        tuple val(sampleID), file(fcontigs), file(depthout)

    output:
    	tuple val(sampleID), file("${sampleID}_bin.*.fa")

	script:
    	"""
    	metabat2 -i $fcontigs -a ${sampleID}_depth.txt -o ${sampleID}_bin -t ${task.cpus} -m ${params.contigsminlength}
    	"""
}

process contigs_to_bins {

	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/Metabat2", mode: 'copy'
	
	input: 
		tuple val(sampleID), file(fafile)

	output:
		tuple val(sampleID), file("${sampleID}_metabat2_contigs_to_bins.tsv"), emit: metabat2_contigs_to_bins

	shell:
    	"""
		grep '>' !{fafile} | tr '>' '\t' | tr -d ':' > !{sampleID}_metabat2_contigs_to_bins.tsv
    	"""
}