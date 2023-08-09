
process KRAKEN2 {

tag "$sampleID"
label 'kraken'
publishDir "${params.outdir}/Kraken/${sampleID}/", mode: 'copy', pattern: "*.txt"

input:
	tuple val(meta), path(reads)

output:
	path(report), emit: krakenreport
	//tuple val(sampleID), file(kraken_log), emit: krakenlog
	tuple val(sampleID), path(report), emit: brackeninput
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	report = sampleID + ".kraken2_report.txt"
	kraken_log = sampleID + "_kraken2.log"

	left_clean = sampleID + "_R1_clean.fastq.gz"
	right_clean = sampleID + "_R2_clean.fastq.gz"
	unpaired_clean = sampleID + "_single_clean.fastq.gz"

    if (!meta.single_end) {  
		"""
		kraken2 --db ${params.kraken2_db} \
			--paired \
			--threads ${task.cpus} \
			--output $kraken_log \
			--report $report ${left_clean} ${right_clean} 
		
		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Kraken2: \$(kraken2 --version | awk 'FNR==1{print \$0}' | sed -e "s/Kraken version //g" )
		END_VERSIONS
	"""
	} else {
		"""
		kraken2 --db ${params.kraken2_db} \
			--threads ${task.cpus} \
			--output $kraken_log \
			--report $report ${unpaired_clean} 
		
		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Kraken2: \$(kraken2 --version | awk 'FNR==1{print \$0}' | sed -e "s/Kraken version //g" )
		END_VERSIONS
		"""	
	}
}

process KRAKEN2MPA {
	label 'default'
	input:
		file(report)

	output:
		path("${report.simpleName}.kraken_mpa.txt"), emit: krakenmpa
		path('versions.yml'), emit: version

	script:
		"""
		kreport2mpa.py -r $report -o ${report.simpleName}.kraken_mpa.txt --percentages --display-header

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
    	END_VERSIONS
		"""
}

process KRAKEN2YAML {
	label 'default'
	input:
		path(reports)

	output:
		path(report_yaml), emit: kraken2yaml
		path('versions.yml'), emit: version

	script:
		report_yaml = "kraken_report_mqc.yaml"
		"""	
		kraken2yaml.pl --outfile $report_yaml

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
    	END_VERSIONS
		"""
}

process KRAKENMERGEREPORTS {
	label 'default'
	publishDir "${params.outdir}/Kraken", mode: 'copy', pattern: "*.txt"

	input:
		path(report)

	output:
		path(report_combined), emit: krakencombined
		path('versions.yml'), emit: version

	script:

		report_combined = "kraken_report_combined.txt"
		"""	
		combine_kreports.py -r ${report.join(" ")} -o $report_combined

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
    	END_VERSIONS
		"""
}

process KRAKENMPAMERGE {
	label 'default'
	publishDir "${params.outdir}/Kraken", mode: 'copy', pattern: "*.txt"

	input:
		path(mpaoutput)

	output:
		path(abundances), emit: krakenmpamerge
		path('versions.yml'), emit: version

	script:
		abundances = "kraken2_mpa_abundances.txt"

		"""
		combine_mpa.py -i ${mpaoutput.join(" ")} -o $abundances

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
    	END_VERSIONS 
		"""
}

process BRACKEN {

	tag "$sampleID"
	label 'bracken'
	publishDir "${params.outdir}/Kraken/${sampleID}/", mode: 'copy', pattern: "*.bracken"

	input:
		tuple val(sampleID), path(report)

	output:
		path(bracken_output), emit: brackenoutput
		path('versions.yml'), emit: version

	script:
		bracken_output = sampleID + ".bracken"
		"""
		bracken -d ${params.kraken2_db} -i ${report} -o ${bracken_output} -r ${params.bracken_length} -l ${params.bracken_level} -t ${params.bracken_threshold}

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version 2>&1 | awk '{print \$2}' )
    	END_VERSIONS 
		"""
}

process BRACKENMERGE {

	label 'bracken'
	publishDir "${params.outdir}/Kraken", mode: 'copy', pattern: "*.txt"

	input:
		path(bracken_output)

	output:
		path(bracken_merged), emit: bracken_merged
		path('versions.yml'), emit: version

	script:
	
	bracken_merged = "bracken_merged.txt"
	"""
	combine_bracken_outputs.py --files ${bracken_output.join(" ")} -o $bracken_merged

	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    Python: \$(python --version 2>&1 | awk '{print \$2}' )
    END_VERSIONS  
	"""
}