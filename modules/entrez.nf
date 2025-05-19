process query_metadata {
	label 'entrez'
	label 'short_run'
	tag "$sampleID"
	scratch params.scratch
	publishDir "${params.outdir}/metadata", mode: 'copy', pattern: "*.csv"

	input:
	tuple val(meta), val(reads)

	output:
	path("${sampleID}_metadata.csv"), optional: true

	script:
	sampleID = meta.id
    
    """
    esearch -db sra -query "${sampleID}" | efetch -format runinfo > ${sampleID}_metadata.csv

    """
}