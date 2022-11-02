process getCountTable {
	label 'default'
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/counttable", mode: 'copy'

	input:
		tuple val(meta), file(finalbam)

	output:
		file("*.txt")

	shell:
		sampleID = meta.id
		"""
		samtools idxstats $finalbam > ${sampleID}_idxstats.txt
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/get_count_table.py ${sampleID}_idxstats.txt > counts_${sampleID}.txt
		"""
}