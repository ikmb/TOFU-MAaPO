process GTDBTK {

	label 'gtdbtk'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/GTDBTK/${sampleID}", mode: 'copy'

	input: 
		tuple val(meta), file(fafile)
		each readygtdbtk

	output:
		path("all.bins.gtdbtk_output/*"), emit: gtdbtk_outputfoldercontent
		tuple val(meta), file("all.bins.gtdbtk_output/gtdbtk.bac120.summary.tsv"), emit: taxonomic_table
		path("versions.yml"),          optional: true, emit: versions
	shell:
		sampleID = meta.id
		"""
		echo "#TRACE n_rows=`tail -n +1 ${fafile} | wc -l`"
		export GTDBTK_DATA_PATH="${params.gtdbtk_reference}"
		gtdbtk classify_wf --cpus ${task.cpus} --genome_dir . --extension fa --out_dir all.bins.gtdbtk_output --pplacer_cpus 1 --skip_ani_screen #/refined_bins

		awk -F "\t" '{ sub(/.*;s__/, "s__", \$2); print \$1 "\t" \$2 }' all.bins.gtdbtk_output/gtdbtk.bac120.summary.tsv > all.bins.gtdbtk_output/parsed_bac120_summary.tsv

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		GTDB-Tk: \$(gtdbtk -version | head -1 | awk '{print \$3}')
		END_VERSIONS
		
		"""
}

    process PREPARE_GTDBTK {

	executor 'local'
    label 'local_run'
    output: 
        val 'true', emit: readystate
	script:

	"""
    if [ ! -d ${params.gtdbtk_reference} ]; then
		mkdir -p ${params.gtdbtk_reference};
    fi
	cd ${params.gtdbtk_reference}

    wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
	tar -xvzf gtdbtk_r207_v2_data.tar.gz
	"""
}