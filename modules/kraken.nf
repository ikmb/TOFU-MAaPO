process KRAKEN2_PE {

tag "$sampleID"
label 'kraken'
publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

input:
tuple val(meta),path(reads),file(unpaired)

output:
path(report), emit: krakenreport
tuple val(sampleID), file(report), emit: brackeninput

script:
sampleID = meta.id
report = sampleID + ".kraken2_report.txt"
kraken_log = sampleID + ".kraken2.log"

"""
kraken2 --db ${params.kraken2_db} --paired --threads ${task.cpus} --output $kraken_log --report $report ${reads[0]} ${reads[1]} 
"""
}

process KRAKEN2_SE {

tag "$sampleID"
label 'kraken'
publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

input:
tuple val(meta),path(reads)

output:
path(report), emit: krakenreport
tuple val(sampleID), file(report), emit: brackeninput

script:
sampleID = meta.id
report = sampleID + ".kraken2_report.txt"
kraken_log = sampleID + ".kraken2.log"

"""
kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --output $kraken_log --report $report ${reads}
"""
}

//output: tuple val(sampleID),file(kraken_log), emit: krakenlog

process KRAKEN2MPA {

input:
file(report)

output:
path("${report.simpleName}.kraken_mpa.txt"), emit: krakenmpa

script:


"""
kreport2mpa.py -r $report -o ${report.simpleName}.kraken_mpa.txt --percentages --display-header
"""
}

process KRAKEN2YAML {

input:
path(reports)

output:
file(report_yaml)

script:

report_yaml = "kraken_report_mqc.yaml"
"""	
kraken2yaml.pl --outfile $report_yaml
"""
}

process KRAKENMERGEREPORTS {
	
publishDir "${params.outdir}/Kraken", mode: 'copy'

input:
path(report)

output:
file(report_combined)

script:

report_combined = "kraken_report_combined.txt"
"""	
combine_kreports.py -r ${report.join(" ")} -o $report_combined
"""
}
// --sample-names ${report.simpleName.join(" ")}
process KRAKENMPAMERGE {

	publishDir "${params.outdir}/Kraken", mode: 'copy'

	input:
	path(mpaoutput)

	output:
	file(abundances)

	script:
	abundances = "kraken2_mpa_abundances.txt"

	"""
		combine_mpa.py -i ${mpaoutput.join(" ")} -o $abundances 
	"""
}

process BRACKEN {

	tag "$sampleID"
	label 'bracken'
	publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

	input:
	tuple val(sampleID), file(report)

	output:
	file(bracken_output)

	script:
	bracken_output = sampleID + ".bracken"
	"""
		bracken -d ${params.kraken2_db} -i ${report} -o ${bracken_output} -r ${params.bracken_length} -l ${params.bracken_level} -t ${params.bracken_threshold}
	"""
}

process BRACKENMERGE {

	publishDir "${params.outdir}/Kraken", mode: 'copy'

	input:
	path(bracken_output)

	output:
	file(bracken_merged)

	script:
	bracken_merged = "bracken_merged.txt"

	"""
		combine_bracken_outputs.py --files ${bracken_output.join(" ")} -o $bracken_merged 
	"""
}