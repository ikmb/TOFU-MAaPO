
process SALMON {
scratch params.scratch
tag "$sampleID"
label 'salmon'
publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "*.sf"
publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "*.log"

input:
	tuple val(meta), path(reads)

output:
	tuple val(sampleID), path(report), path(salmon_log), emit: salmonout
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	report = sampleID + ".quant.sf"
	salmon_log = sampleID + "_salmon.log"

    if (!meta.single_end) {
	"""
    salmon quant -i ${params.salmon_db} \
        -l IU \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
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
	} else {
	"""
	salmon quant -i ${params.salmon_db} \
        -l IU \
        -1 ${reads[0]} \
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
}

