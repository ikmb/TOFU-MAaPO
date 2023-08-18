process getCountTable {
	label 'default'
	cache 'lenient'
	tag "$sampleID"
	publishDir "${params.outdir}/counttable/${sampleID}", mode: 'copy'

	input:
		tuple val(meta), path(finalbam), path(mappingbam_index)
	
	output:
		path("*.txt"),							emit: counttableoutput
		path("versions.yml"),	optional: true, emit: versions

	shell:
		sampleID = meta.id
		"""
		samtools idxstats $finalbam > ${sampleID}_idxstats.txt
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/get_count_table.py ${sampleID}_idxstats.txt > counts_${sampleID}.txt

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
      	Python: \$(/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 --version | sed -e "s/Python //g" )
		samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
		END_VERSIONS
		
		"""
}