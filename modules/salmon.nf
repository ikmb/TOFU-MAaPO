
process KRAKEN2 {

tag "$sampleID"
label 'kraken'
publishDir "${params.outdir}/Kraken/${sampleID}/", mode: 'copy'

input:
	tuple val(meta), path(reads)

output:
	path(report), emit: krakenreport
	//tuple val(sampleID), file(kraken_log), emit: krakenlog
	tuple val(sampleID), path(report), emit: brackeninput
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	report = sampleID + ".quant.sf"
	salmon_log = sampleID + "_salmon.log"

    if (!params.single_end) {  
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
    salmon: \$(salmon --version 2>&1 | | sed -e "s/salmon //g" )
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
    Kraken2: \$(kraken2 --version | awk 'FNR==1{print \$0}' | sed -e "s/Kraken version //g" )
    END_VERSIONS
	"""	
	}
}

