process GTDBTK {

	label 'gtdbtk'
	//scratch params.scratch
	scratch false
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/gtdbtk", mode: 'copy'

	input: 
	    tuple val(meta), file(fafile)

	output:
	    file("all.bins.gtdbtk_output/*")
	
	shell:
		sampleID = meta.id
	    """
	    export GTDBTK_DATA_PATH=${params.gtdbtk_reference}
	    gtdbtk classify_wf --cpus ${task.cpus} --genome_dir ./refined_bins --extension fa --out_dir all.bins.gtdbtk_output --pplacer_cpus 1

	    gawk -F "\t" '{ sub(/.*;s__/, "s__", \$2); print \$1 "\t" \$2 }' all.bins.gtdbtk_output/gtdbtk.bac120.summary.tsv > all.bins.gtdbtk_output/parsed_bac120_summary.tsv
	    """
}