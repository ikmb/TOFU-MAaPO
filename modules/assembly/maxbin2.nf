process MAXBIN2 {

	label 'maxbin2'
	//scratch params.scratch
    scratch false
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/maxbin2", mode: 'copy'

	input:
	    tuple val(sampleID), file(fcontigs), file(depthout)
	output:
		//file("*"), emit: all
		tuple val(sampleID), file(maxbin2_contigs_to_bin), emit: contigs_to_bin

	script:
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
