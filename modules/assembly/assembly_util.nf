process getCountTable {
	
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/counttable", mode: 'copy'

	input:
		tuple val(sampleID), file(finalbam)

	output:
		file("*.txt")

	shell:
		"""
		samtools idxstats $finalbam > ${sampleID}_idxstats.txt
		python ${baseDir}/bin/get_count_table.py ${sampleID}_idxstats.txt > counts_${sampleID}.txt
		"""
}