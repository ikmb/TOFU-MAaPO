process BINCOVERAGE_PERSAMPLE{
	label 'default'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/MAG_abundance/${sampleID}", mode: 'copy', pattern: "*.tbl"

	input:
		tuple val(meta), file(depthout), file(refined_contigs_to_bins), file(taxonomic_table)
	output:
		path(output_file), emit: abundancetable
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id
		output_file = sampleID + '_abundance_table.tbl'
		"""
		Rscript ${baseDir}/bin/coverage_table.R $depthout 150 $refined_contigs_to_bins $taxonomic_table $sampleID $output_file
		#echo "suppressMessages(library(tidyverse));suppressMessages(library(data.table));sessionInfo()" | R --slave > R_version.txt

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		R: \$(Rscript --version 2>&1 | awk '{print \$5}')
		END_VERSIONS
		"""
}

process MERGE_MAG_ABUNDANCE{
	label 'default'
	scratch params.scratch
	publishDir "${params.outdir}/MAG_abundance", mode: 'copy', pattern: "*.tbl"
	publishDir "${params.outdir}/MAG_abundance", mode: 'copy', pattern: "*.png"

	input:
		file(input_files)
	output:
		file(output_file)
		file(output_plot)
		path("versions.yml"),          optional: true, emit: versions

	script:
		output_file = 'merged_MAG_tpm_abundance.tbl'
		output_plot = 'phylum_rel_abudance_plot.png'
		"""
		Rscript ${baseDir}/bin/merge_MAG_abundance.R

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		R: \$(Rscript --version 2>&1 | awk '{print \$5}')
		END_VERSIONS
		"""
}