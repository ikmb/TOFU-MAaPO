process MAXBIN2 {

	label 'maxbin2'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/maxbin2/${sampleID}", mode: 'copy', enabled: params.publish_rawbins
	
	input:
	    tuple val(meta), file(fcontigs), file(depthout)
	output:
		tuple val(meta), file(maxbin2_contigs_to_bin), emit: contigs_to_bin

	script:
		sampleID = meta.id
    	abundance_table = sampleID + '.abu'
    	maxbin2output = 'maxbin2_out/' + sampleID + '.maxbin2'
		maxbin2_contigs_to_bin = sampleID + '_maxbin2_contigs_to_bin.tsv'

	"""
        awk '{if(NR>1) print \$1"\\t"\$3}' $depthout > $abundance_table
        
		mkdir -p maxbin2_out

		run_MaxBin.pl \
            -contig $fcontigs \
            -out $maxbin2output \
            -abund $abundance_table \
            -thread ${task.cpus}	 \
            -min_contig_length 2000
		
		grep '>' maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\\ /' > $maxbin2_contigs_to_bin
	"""
}
