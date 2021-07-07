process KRAKEN2 {

label 'kraken'

publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

input:
tuple val(sampleID),file(left),file(right),file(unpaired)

output:
path(report), emit: krakenreport

script:
report = sampleID + ".kraken2_report.txt"
kraken_log = sampleID + ".kraken2.log"

"""
kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --output $kraken_log --report $report $left $right
"""
}
//output: tuple val(sampleID),file(kraken_log), emit: krakenlog

process KR2MPA {

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