process FORMATTING_CONTIG_TO_BIN {
	label 'default'
	scratch params.scratch
	tag "${sampleID}_${meta.assembler}"

	input:
		tuple val(meta), file(contigs_bin_tables)
	output:
		tuple val(meta), file(merged_contigs_to_bin), emit: formatted_contigs_to_bin

	script:
		sampleID = meta.id

		merged_contigs_to_bin = sampleID + '_magscot_contigs_to_bin_merged.tsv'
		"""
		cat ${contigs_bin_tables.join(" ")} > $merged_contigs_to_bin
		"""
	stub:
        sampleID = meta.id
		merged_contigs_to_bin = sampleID + '_magscot_contigs_to_bin_merged.tsv'
		"""
		touch $merged_contigs_to_bin

		echo "FORMATTING_CONTIG_TO_BIN_stub" > versions.yml
		"""
}

process MARKER_IDENT {
	label 'magscot'
	scratch params.scratch
	tag "${sampleID}_${assembler}"
	publishDir "${params.outdir}/magscot/${sampleID}_${assembler}", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), file(depthout), file(formatted_contigs_to_bin)
	output:
		tuple val(meta), file(samplehmm), emit: hmm_output
		path("versions.yml"),          optional: true, emit: versions
	script:
		sampleID = meta.id
		assembler = meta.assembler
		samplehmm = sampleID + '.hmm'
		samplepfam = sampleID + '.pfam'
		sampleptigr = sampleID + '.tigr'
		sampleprodigalfaa = sampleID + '.prodigal.faa'
		sampleprodigalffn = sampleID + '.prodigal.ffn'
		sample_tmp = sampleID + '_tmp'
		"""
		### ORF detection with prodigal
		cat $fcontigs | prodigal -p meta -a $sampleprodigalfaa -d $sampleprodigalffn -o $sample_tmp

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
	stub:
		sampleID = meta.id
		assembler = meta.assembler
		samplehmm = sampleID + '.hmm'
		"""
		touch $samplehmm

		echo "MARKER_IDENT_stub" > versions.yml
		"""
}

process MAGSCOT {

	label 'magscot'
	scratch params.scratch
	errorStrategy  { task.exitStatus in [42] ? 'ignore' : task.attempt <= maxRetries  ? 'retry' : 'ignore' }
	tag "${sampleID}_${meta.assembler}"
	publishDir "${params.outdir}/magscot/${sampleID}_${meta.assembler}", mode: 'copy'

	input:
		tuple val(coassemblygroup), val(meta), file(formatted_contigs_to_bin), file(samplehmm), file(fcontigs_filtered)
	output:
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered), emit: refined_contigs_to_bins
		tuple val(meta), file(refined_contigs_to_bins), emit: contigs_to_bins_table
		tuple val(meta), file(stats_outfile), emit: stats_outfile_table
		tuple val(meta), file(full_stats), emit: full_bin_stats
		path("versions.yml"),          optional: true, emit: versions
	script:
		sampleID = meta.id
		outputname = meta.id + '_' + meta.assembler
		refined_contigs_to_bins = outputname + '.refined.contig_to_bin.out'
		stats_outfile = outputname + '.refined.out'
		full_stats = outputname + '.scores.out'
	"""
		#If the input couldn't be binned with any binner, we ignore that specific error and create a specific error code to ignore it with nextflow
		set +e
		{
		Rscript /opt/MAGScoT.R \
			-i $formatted_contigs_to_bin \
			--hmm $samplehmm \
			-o $outputname \
			-s ${params.magscot_min_sharing} \
			-a ${params.magscot_score_a} \
			-b ${params.magscot_score_b} \
			-c ${params.magscot_score_c} \
			-t ${params.magscot_threshold} \
			-m ${params.magscot_min_markers} \
			-n ${params.magscot_iterations}  > log.txt

		if grep -q "HMM input file should have three tab-separated colums without a header: Protein ID, Marker ID, e-value" log.txt; then
            # Set the error code to 42
            exit_code=42
        else
            # If no error or different error, continue with the script
            echo "No specific error encountered, continuing with the script."
            exit_code=\$?
        fi

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		MAGScoT: 1.0.0
		R: \$(Rscript --version | awk '{print \$4}')
		END_VERSIONS

		}
		set -e

		exit \$exit_code
	"""
	stub:
		sampleID = meta.id
		outputname = meta.id + '_' + meta.assembler
		refined_contigs_to_bins = outputname + '.refined.contig_to_bin.out'
		stats_outfile = outputname + '.refined.out'
		full_stats = outputname + '.scores.out'
		"""
		touch $refined_contigs_to_bins
		touch $stats_outfile
		touch $full_stats

		echo "MAGSCOT_stub" > versions.yml
		"""
}

process EXTRACT_REFINED_BINS {
	label 'default'
	scratch params.scratch
	tag "${sampleID}_${meta.assembler}"
	publishDir "${params.outdir}/magscot/${sampleID}_${meta.assembler}", mode: 'copy'

	input:
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered)
	output:
		tuple val(meta), file("refined_bins/*"), emit: refined_bins
		path("versions.yml"),          optional: true, emit: versions
		//tuple val(meta), file("refined_bins"), emit: refined_bins_folder
	script:
		sampleID = meta.id
	"""
		mkdir -p refined_bins
		
		cat $refined_contigs_to_bins | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print \$2","\$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $fcontigs_filtered /dev/stdin  --output_path refined_bins

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS
		
	"""
	stub:
		sampleID = meta.id
		"""
		mkdir refined_bins
		touch refined_bins/${sampleID}_stub.fa

		echo "EXTRACT_REFINED_BINS_stub" > versions.yml
		"""
}
