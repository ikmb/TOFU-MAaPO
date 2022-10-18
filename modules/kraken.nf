
process KRAKEN2 {

tag "$sampleID"
label 'kraken'
publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

input:
	tuple val(meta), path(reads)

output:
	path(report), emit: krakenreport
	//tuple val(sampleID), file(kraken_log), emit: krakenlog
	tuple val(sampleID), file(report), emit: brackeninput

script:
	sampleID = meta.id
	report = sampleID + ".kraken2_report.txt"
	kraken_log = sampleID + "_kraken2.log"

	left_clean = sampleID + "_R1_clean.fastq.gz"
	right_clean = sampleID + "_R2_clean.fastq.gz"
	unpaired_clean = sampleID + "_single_clean.fastq.gz"

    if (!params.single_end) {  
	"""
	kraken2 --db ${params.kraken2_db} \
		--paired \
		--threads ${task.cpus} \
		--output $kraken_log \
		--report $report ${left_clean} ${right_clean} 
	"""
	} else {
	"""
	kraken2 --db ${params.kraken2_db} \
		--threads ${task.cpus} \
		--output $kraken_log \
		--report $report ${unpaired_clean} 
	"""	
	}
}

process KRAKEN2MPA {
	label 'default'
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
	label 'default'
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
	label 'default'
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
	label 'default'
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

	label 'bracken'
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