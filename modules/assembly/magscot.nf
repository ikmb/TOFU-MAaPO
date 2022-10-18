process FORMATTING_CONTIG_TO_BIN {
	label 'default'
	//scratch params.scratch
    scratch false
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/vamb", mode: 'copy'

	input:
        tuple val(meta), file(vamb_cluster_table), file(metabat2_cluster_table), file(maxbin2_cluster_table), file(concoct_cluster_table)

	output:
		tuple val(meta), file(formatted_contigs_to_bin), emit: formatted_contigs_to_bin

	script:
		sampleID = meta.id
		formatted_contigs_to_bin = sampleID + '_contigs_to_bin.tsv'
	"""
    	gawk '{print \$1"\t"\$2"\tvamb"}'  $vamb_cluster_table > $formatted_contigs_to_bin
		gawk '{print \$1"\t"\$2"\tmetabat2"}'  $metabat2_cluster_table >> $formatted_contigs_to_bin
		gawk '{print \$1"\t"\$2"\tmaxbin2"}'  $maxbin2_cluster_table >> $formatted_contigs_to_bin
		gawk '{print \$1"\t"\$2"\tconcoct"}'  $concoct_cluster_table >> $formatted_contigs_to_bin

	"""
}

process MARKER_IDENT {

	label 'magscot'
	//scratch params.scratch
    scratch false
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/magscot", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), file(depthout), file(formatted_contigs_to_bin)
	output:
		//tuple val(sampleID), file(formatted_contigs_to_bin), emit: formatted_contigs_to_bin
		tuple val(meta), file(samplehmm), emit: hmm_output
	script:
		sampleID = meta.id

		samplehmm = sampleID + '.hmm'
		samplepfam = sampleID + '.pfam'
		sampleptigr = sampleID + '.tigr'
		sampleprodigalfaa = sampleID + '.prodigal.faa'
		sampleprodigalffn = sampleID + '.prodigal.ffn'
		sample_tmp = sampleID + '_tmp'
	"""
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

		cat ${sampleID}.hmm.tigr.hit.out | grep -v "^#" | awk '{print \$1"\t"\$3"\t"\$5}' > $sampleptigr
		cat ${sampleID}.hmm.pfam.hit.out | grep -v "^#" | awk '{print \$1"\t"\$4"\t"\$5}' > $samplepfam
		cat $samplepfam $sampleptigr > $samplehmm
	"""
}

process MAGSCOT {

	label 'magscot'
	scratch params.scratch
    //scratch false
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/magscot", mode: 'copy'

	input:
		tuple val(meta), file(formatted_contigs_to_bin), file(samplehmm), file(fcontigs_filtered)
	output:
		//tuple val(sampleID), file(formatted_contigs_to_bin), emit: formatted_contigs_to_bin
		file("*")
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered), emit: refined_contigs_to_bins
	script:
		sampleID = meta.id
		refined_contigs_to_bins = sampleID + '.refined.contig_to_bin.out'
	"""
		Rscript /opt/MAGScoT.R -i $formatted_contigs_to_bin --hmm $samplehmm -o $sampleID

	"""
}

process EXTRACT_REFINED_BINS {
	label 'default'
	//scratch params.scratch
    scratch false
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/magscot", mode: 'copy'

	input:
		tuple val(meta), file(refined_contigs_to_bins), file(fcontigs_filtered)
	output:
		tuple val(meta), file("refined_bins"), emit: refined_bins
	script:

	"""
		mkdir -p refined_bins
		
		cat $refined_contigs_to_bins | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print \$2","\$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $fcontigs_filtered /dev/stdin  --output_path refined_bins
	"""
}