process MAXBIN2 {

	label 'maxbin2'
	scratch params.scratch
	tag "${sampleID}_${meta.assembler}"
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
		awk '{if(NR>1) print \$1"\\t"\$3}' $depthout > $abundance_table
		
		mkdir -p maxbin2_out

		#If the input cannot be assembled, we ignore that specific error and create an empty output
		set +e
		{
		output=\$(run_MaxBin.pl \
			-contig $fcontigs \
			-out $maxbin2output \
			-abund $abundance_table \
			-thread ${task.cpus}	 \
			-min_contig_length ${params.contigsminlength} \
			-preserve_intermediate )
		}
        set -e

		if [ \$? -ne 0 ] && echo "\$output" | grep -q "Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop."; then
            # Set the error code back to 0
            exit_code=0
        else
            # If no error or different error, continue with the script
            echo "Not encountered the specific error, continuing with the script."
            exit_code=\$?
        fi

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
		
		# Exit with the modified exit code
		exit \$exit_code
	"""
	stub:
		sampleID = meta.id
		maxbin2_contigs_to_bin = sampleID + '_maxbin2_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_maxbin2_magscot_contigs_to_bin.tsv'
		"""
		touch $maxbin2_contigs_to_bin
		touch $formatted_contigs_to_bin

		echo "MAXBIN2_stub" > versions.yml
		"""
}
