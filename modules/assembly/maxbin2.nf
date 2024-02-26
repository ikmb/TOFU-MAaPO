process MAXBIN2 {

	label 'maxbin2'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/maxbin2/${sampleID}", mode: 'copy', enabled: params.publish_rawbins
	
	input:
		tuple val(meta), file(fcontigs), file(depthout)
	output:
		tuple val(meta), file(maxbin2_contigs_to_bin), optional: true, emit: contigs_to_bin
		tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
		path("versions.yml"),          optional: true, emit: versions
	script:
		sampleID = meta.id
		abundance_table = sampleID + '.abu'
		maxbin2output = 'maxbin2_out/' + sampleID + '.maxbin2'
		maxbin2_contigs_to_bin = sampleID + '_maxbin2_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_maxbin2_magscot_contigs_to_bin.tsv'
	"""
		echo "#TRACE n_rows=`tail -n +1 ${fcontigs} | wc -l`"

		awk '{if(NR>1) print \$1"\\t"\$3}' $depthout > $abundance_table
		
		mkdir -p maxbin2_out

		run_MaxBin.pl \
			-contig $fcontigs \
			-out $maxbin2output \
			-abund $abundance_table \
			-thread ${task.cpus}	 \
			-min_contig_length ${params.contigsminlength} \
			-preserve_intermediate

		rm maxbin2_out/${sampleID}.maxbin2.contig.tmp.*

		set +e
		{
			grep '>' maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\\ /' > $maxbin2_contigs_to_bin

			awk '{print \$1"\t"\$2"\tmaxbin2"}' $maxbin2_contigs_to_bin > $formatted_contigs_to_bin
		}
		set -e

		cat <<-END_VERSIONS> versions.yml
		"${task.process}":
		maxbin2: \$(run_MaxBin.pl -v 2>&1 | tail -1 | sed -e 's/MaxBin //g')
		END_VERSIONS
		
	"""
}
