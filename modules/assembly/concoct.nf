process CONCOCT {
	label 'concoct'
    label 'long_run'
    cache 'lenient'
    errorStrategy { task.exitStatus == 43 ? 'ignore' : task.attempt <= 2 ? 'retry' : 'finish' } //this would not publish failed logs
	scratch params.scratch
	tag "$sampleID"
    publishDir "${params.outdir}/concoct", mode: 'copy', enabled: params.publish_rawbins,
        saveAs: { filename -> "${meta.id}/${filename}" }

    publishDir "${params.outdir}/concoct",
        mode: 'copy',
        pattern: "*_concoct_failed.txt",
        saveAs: { filename -> "${meta.id}/${filename}" }

	input:
        tuple val(meta), file(fcontigs), file(depthout), file(mappingbam), file(mappingbam_index)
	
    output:
		tuple val(meta), file(concoct_contigs_to_bin), optional: true, emit: contigs_to_bin
        tuple val(meta), file(formatted_contigs_to_bin), optional: true, emit: magscot_contigbinlist
		tuple val(meta), file("${sampleID}_concoct_failed.txt"), optional: true, emit: concoct_failed
        path("versions.yml"),          optional: true, emit: versions
    script:
        sampleID = meta.id
        bed_file = sampleID + '.bed'
        concoct_contigs_to_bin = sampleID + '_concoct_contigs_to_bin.tsv'
        formatted_contigs_to_bin = sampleID + '_concoct_magscot_contigs_to_bin.tsv'
        failed_log = sampleID + '_concoct_failed.txt'
        """
        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/cut_up_fasta.py $fcontigs -c 10000 -o 0 --merge_last -b $bed_file > ${sampleID}.filtered.10k.fna

        samtools index $mappingbam

        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/concoct_coverage_table.py $bed_file $mappingbam > ${sampleID}.coverage_table.tsv
        
        if concoct \
            --composition_file ${sampleID}.filtered.10k.fna \
            --coverage_file ${sampleID}.coverage_table.tsv \
            -c 1000 \
            -r 151 \
            -t ${task.cpus} \
            -l ${params.contigsminlength} \
            -s 1234 \
            -i 500 \
            -b ${sampleID}
        then
            :
        else
            concoct_exit=\$?

            echo "CONCOCT exited with: \$concoct_exit" >&2
            echo "Checking CONCOCT log: ${sampleID}_log.txt" >&2

            if grep -Fq "Not enough contigs pass the threshold filter" ${sampleID}_log.txt; then
                {
                    echo "WARNING: CONCOCT skipped sample ${sampleID}: not enough contigs pass the threshold filter."
                    echo
                    cat ${sampleID}_log.txt
                } > ${sampleID}_concoct_failed.txt

                exit 0
            else
                echo "CONCOCT encountered an error" >&2
                exit "\$concoct_exit"
            fi


        fi

        /opt/conda/envs/ikmb-metagenome-1.2/bin/python3.10 /opt/conda/envs/ikmb-metagenome-1.2/bin/merge_cutup_clustering.py ${sampleID}_clustering_gt2000.csv > ${sampleID}_clustering_merged.csv

        awk -F',' -v sample=${sampleID} '{if(NR>1) print sample"_concoct_bin_"\$2".fasta\t"\$1}'  ${sampleID}_clustering_merged.csv > $concoct_contigs_to_bin

        awk '{print \$1"\t"\$2"\tconcoct"}'  $concoct_contigs_to_bin > $formatted_contigs_to_bin

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		concoct: \$(concoct -v |  awk '{print \$2}')
		END_VERSIONS
        """
    stub:
        sampleID = meta.id
        bed_file = sampleID + '.bed'
        concoct_contigs_to_bin = sampleID + '_concoct_contigs_to_bin.tsv'
        formatted_contigs_to_bin = sampleID + '_concoct_magscot_contigs_to_bin.tsv'

		"""
		touch $concoct_contigs_to_bin
        touch $formatted_contigs_to_bin

		echo "CONCOCT_stub" > versions.yml
		"""
}