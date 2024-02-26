process FORMATTING_CONTIG_TO_BIN {
	label 'default'
	scratch params.scratch
	tag "$sampleID"

	input:
		tuple val(meta), file(contigs_bin_tables)
	output:
		tuple val(meta), file(merged_contigs_to_bin), emit: formatted_contigs_to_bin

	script:
		sampleID = meta.id

		merged_contigs_to_bin = sampleID + '_magscot_contigs_to_bin_merged.tsv'
		"""
		echo "#TRACE n_rows=`tail -n +1 ${contigs_bin_tables} | wc -l`"
		cat ${contigs_bin_tables.join(" ")} > $merged_contigs_to_bin
		"""
}

process MARKER_IDENT {

	label 'magscot'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/magscot/${sampleID}", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), file(depthout), file(formatted_contigs_to_bin)
	output:
		tuple val(meta), file(samplehmm), emit: hmm_output
		path("versions.yml"),          optional: true, emit: versions
	script:
		sampleID = meta.id

		samplehmm = sampleID + '.hmm'
		samplepfam = sampleID + '.pfam'
		sampleptigr = sampleID + '.tigr'
		sampleprodigalfaa = sampleID + '.prodigal.faa'
		sampleprodigalffn = sampleID + '.prodigal.ffn'
		sample_tmp = sampleID + '_tmp'
		"""
		echo "#TRACE n_rows=`tail -n +1 ${fcontigs} | wc -l`"
		### ORF detection with prodigal
		cat $fcontigs | prodigal -p meta -a $sampleprodigalfaa -d $sampleprodigalffn -o $sample_tmp

		### ALTERNATIVE: Fast parallel ORF detection with prodigal
		# mkdir -p tmp_workfolder
		# zcat example.contigs.fasta.gz | parallel -j 8 --block 999k --recstart '>' --pipe prodigal -p meta -a tmp_workfolder/example.{#}.faa -d # tmp_workfolder/example.{#}.ffn -o tmpfile
		# cat tmp_workfolder/example.*.faa > example.prodigal.faa
		# cat tmp_workfolder/example.*.ffn > example.prodigal.ffn
		# rm -r tmp_workfolder tmpfile

		### annotation of protein sequences using HMMer and GTDBtk r207 marker genes
		hmmsearch \
			-o ${sampleID}.hmm.tigr.out \
			--tblout ${sampleID}.hmm.tigr.hit.out \
			--noali \
			--notextw \
			--cut_nc \
			--cpu ${task.cpus} \
			/opt/hmm/gtdbtk_rel207_tigrfam.hmm \
			$sampleprodigalfaa
		hmmsearch \
			-o ${sampleID}.hmm.pfam.out \
			--tblout ${sampleID}.hmm.pfam.hit.out \
			--noali \
			--notextw \
			--cut_nc \
			--cpu ${task.cpus} \
			/opt/hmm/gtdbtk_rel207_Pfam-A.hmm \
			$sampleprodigalfaa
		set +e
		{
		cat ${sampleID}.hmm.tigr.hit.out | grep -v "^#" | awk '{print \$1"\t"\$3"\t"\$5}' > $sampleptigr
		cat ${sampleID}.hmm.pfam.hit.out | grep -v "^#" | awk '{print \$1"\t"\$4"\t"\$5}' > $samplepfam

		cat $samplepfam $sampleptigr > $samplehmm
		}
		set -e

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		hhmsearch: \$(hmmsearch -h 2>&1 | awk 'NR==2{print \$3}')
		Prodigal: \$(prodigal -v 2>&1 | awk 'NR==2{print\$2}' | sed -e 's/V//g' | sed -e 's/://g')
		R: \$(Rscript --version | awk '{print \$4}')
		END_VERSIONS

		"""
}

process MAGSCOT {

	label 'magscot'
	scratch params.scratch
	errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
	tag "$sampleID"
	publishDir "${params.outdir}/magscot/${sampleID}", mode: 'copy'

	input:
		tuple val(coassemblygroup), val(meta), file(formatted_contigs_to_bin), file(samplehmm), file(fcontigs_filtered)
	output:
		//file("*"), emit: all_files
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered), emit: refined_contigs_to_bins
		tuple val(meta), file(refined_contigs_to_bins), emit: contigs_to_bins_table
		tuple val(meta), file(stats_outfile), emit: stats_outfile_table
		path("versions.yml"),          optional: true, emit: versions
	script:
		sampleID = meta.id
		refined_contigs_to_bins = sampleID + '.refined.contig_to_bin.out'
		stats_outfile = sampleID + '.refined.out'
		full_stats = sampleID + '.scores.out'
	"""
		echo "#TRACE n_rows=`tail -n +1 ${formatted_contigs_to_bin} | wc -l`"
		Rscript /opt/MAGScoT.R -i $formatted_contigs_to_bin --hmm $samplehmm -o $sampleID -s ${params.magscot_min_sharing}

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		MAGScoT: 1.0.0
		R: \$(Rscript --version | awk '{print \$4}')
		END_VERSIONS

	"""
}

process EXTRACT_REFINED_BINS {
	label 'default'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/magscot/${sampleID}", mode: 'copy'

	input:
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered)
	output:
		tuple val(meta), file("refined_bins/*"), emit: refined_bins
		path("versions.yml"),          optional: true, emit: versions
		//tuple val(meta), file("refined_bins"), emit: refined_bins_folder
	script:
		sampleID = meta.id
	"""
		echo "#TRACE n_rows=`tail -n +1 ${refined_contigs_to_bins} | wc -l`"
		mkdir -p refined_bins
		
		cat $refined_contigs_to_bins | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print \$2","\$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $fcontigs_filtered /dev/stdin  --output_path refined_bins

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS
		
	"""
}
