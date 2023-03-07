process METABAT {

    label 'metabat2'
    scratch params.scratch
    tag "$sampleID"
    publishDir "${params.outdir}/${sampleID}/Metabat2", mode: 'copy', enabled: params.publish_rawbins

    input: 
        tuple val(meta), file(fcontigs), file(depthout)

    output:
    	tuple val(meta), file("${sampleID}_bin.*.fa")

	script:
		sampleID = meta.id
    	"""
    	metabat2 -i $fcontigs -a ${sampleID}_depth.txt -o ${sampleID}_bin -t ${task.cpus} -m ${params.contigsminlength}
    	"""
}

process contigs_to_bins {

	label 'default'
	scratch params.scratch
	tag "$sampleID"
	publishDir {"${params.outdir}/${sampleID}/Metabat2"}, mode: 'copy', enabled: params.publish_rawbins as boolean
	
	input: 
		tuple val(meta), file(fafile)

	output:
		tuple val(meta), file("${sampleID}_metabat2_contigs_to_bins.tsv"), emit: metabat2_contigs_to_bins

	script:
		sampleID = meta.id
		
		if (fafile instanceof List && fafile.size() > 1) {
        	"""
        	echo "Variable contains multiple files"
			grep '>' ${fafile} | tr '>' '\t' | tr -d ':' > ${sampleID}_metabat2_contigs_to_bins.tsv
       		"""
    	} else {
			"""
			echo "Variable contains a single file"
			grep '>' ${fafile} | tr '>' '\t' | tr -d ':' | awk -F'\t' '{OFS="\t"; if (\$1=="") \$1="${fafile}"; print \$0 }' > ${sampleID}_metabat2_contigs_to_bins.tsv
    		"""
    	}
    	
}