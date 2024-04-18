
process SALMON {
scratch params.scratch
tag "$sampleID"
label 'salmon'
publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "*.sf"
publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "*.log"

input:
	tuple val(meta), path(reads)

output:
	tuple val(sampleID), path(report), path(salmon_log), emit: salmonallout
	path(report), emit: salmonout
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	report = sampleID + ".quant.sf"
	salmon_log = sampleID + "_salmon.log"
	def input = meta.single_end ? "-1 ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}" 

	"""
	salmon quant -i ${params.salmon_db} \
		-l IU \
		$input \
		--validateMappings \
		-o . \
		-p ${task.cpus} \
		--meta
	
	mv quant.sf $report
	mv logs/salmon_quant.log $salmon_log

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	salmon: \$(salmon --version 2>&1 | sed -e "s/salmon //g" )
	END_VERSIONS

	"""
}

process SALMON_merge {
	label 'default'
	publishDir "${params.outdir}/SALMON", mode: 'copy', pattern: "*.tbl"

	input:
		path(report)

	output:
		path(abundances),          optional: true, emit: salmonmerge
		path('versions.yml'), emit: version

	script:
		def args = params.salmon_reference ? "${params.salmon_reference}" : ""
		abundances = "salmon_merged_TPM.tbl"
		"""	
		Rscript ${baseDir}/bin/salmon_merging.R ${args}

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		R: \$(Rscript --version 2>&1 | awk '{print \$5}')
		END_VERSIONS
		"""	
}

