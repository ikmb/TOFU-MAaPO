process CONCOCT {
	label 'concoct'
    cache 'lenient'
	scratch params.scratch
	tag "$sampleID"
	publishDir {"${params.outdir}/concoct/${sampleID}"}, mode: 'copy', enabled: params.publish_rawbins

	input:
        tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)
	
    output:
		tuple val(meta), file(concoct_contigs_to_bin), optional: true, emit: contigs_to_bin
        tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
        path("versions.yml"),          optional: true, emit: versions
    script:
        sampleID = meta.id
        bed_file = sampleID + '.bed'
        concoct_contigs_to_bin = sampleID + '_concoct_contigs_to_bin.tsv'
        formatted_contigs_to_bin = sampleID + 'concoct_magscot_contigs_to_bin.tsv'
        """
        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/cut_up_fasta.py $fcontigs -c 10000 -o 0 --merge_last -b $bed_file > ${sampleID}.filtered.10k.fna

        samtools index $mappingbam

        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/concoct_coverage_table.py $bed_file $mappingbam > ${sampleID}.coverage_table.tsv

        concoct \
            --composition_file ${sampleID}.filtered.10k.fna \
            --coverage_file ${sampleID}.coverage_table.tsv \
            -c 1000 \
            -r 151 \
            -t ${task.cpus} \
            -l ${params.contigsminlength} \
            -s 1234 \
            -i 500 \
            -b ${sampleID}

        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/merge_cutup_clustering.py ${sampleID}_clustering_gt2000.csv > ${sampleID}_clustering_merged.csv

        awk -F',' -v sample=${sampleID} '{if(NR>1) print sample"_concoct_bin_"\$2".fasta\t"\$1}'  ${sampleID}_clustering_merged.csv > $concoct_contigs_to_bin

        awk '{print \$1"\t"\$2"\tconcoct"}'  $concoct_contigs_to_bin > $formatted_contigs_to_bin

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		concoct: \$(concoct -v |  awk '{print \$2}')
		END_VERSIONS
        """
}