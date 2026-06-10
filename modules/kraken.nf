
process KRAKEN2 {

tag "${meta.id}"
label 'kraken'
label 'long_run'
publishDir "${params.outdir}/Kraken", mode: 'copy', pattern: "*.txt",
        saveAs: { filename -> "${meta.id}/${filename}" }


input:
	tuple val(meta), path(reads)

output:
	path("${meta.id}.kraken2_report.txt"), emit: krakenreport
	tuple val(meta.id), path("${meta.id}.kraken2_report.txt"), emit: brackeninput
	path('versions.yml'), emit: version

script:
	def sampleID = meta.id
	def report = sampleID + ".kraken2_report.txt"
	def kraken_log = sampleID + "_kraken2.log"

	def left_clean = sampleID + "_R1_clean.fastq.gz"
	def right_clean = sampleID + "_R2_clean.fastq.gz"
	def unpaired_clean = sampleID + "_single_clean.fastq.gz"

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
	stub:
	def sampleID = meta.id
	def report = sampleID + ".kraken2_report.txt"
	def kraken_log = sampleID + "_kraken2.log"
	"""
	cat > ${report} <<EOF
# Kraken2 stub report for ${sampleID}
100.00\t1\t1\tU\t0\tunclassified
0.00\t0\t0\tR\t1\troot
EOF

cat > ${kraken_log} <<EOF
[stub] kraken2 executed for ${sampleID}
EOF

cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Kraken2: stub
END_VERSIONS
"""
}

process KRAKEN2MPA {
	label 'default'
	label 'short_run'
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
		stub:
		"""
		cat > ${report.simpleName}.kraken_mpa.txt <<EOF
#Classification\tstub_sample
k__Bacteria\t100.0
EOF

		cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
		"""
}

process KRAKEN2YAML {
	label 'default'
	label 'short_run'
	input:
		path(reports)

	output:
		path("kraken_report_mqc.yaml"), emit: kraken2yaml
		path('versions.yml'), emit: version

	script:
		def report_yaml = "kraken_report_mqc.yaml"
		"""	
		kraken2yaml.pl --outfile $report_yaml

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
    	END_VERSIONS

		"""
	stub:
		def report_yaml = "kraken_report_mqc.yaml"
		"""
		cat > ${report_yaml} <<EOF
id: kraken2
section_name: Kraken2
description: Stub Kraken2 summary
plot_type: table
pconfig:
  id: kraken2_stub_table
data:
  stub_sample:
    classified: 1
    unclassified: 0
EOF

		cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
		"""
}

process KRAKENMERGEREPORTS {
	label 'default'
	label 'short_run'
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
	stub:
	report_combined = "kraken_report_combined.txt"
	"""
	cat > ${report_combined} <<EOF
# Kraken2 combined stub report
100.00\t1\t1\tU\t0\tunclassified
0.00\t0\t0\tR\t1\troot
EOF

	cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
	"""
}

process KRAKENMPAMERGE {
	label 'default'
	label 'short_run'
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
	stub:
		abundances = "kraken2_mpa_abundances.txt"
		"""
		cat > ${abundances} <<EOF
Classification\tstub_sample
k__Bacteria\t100.0
EOF

		cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
		"""
}

process BRACKEN {

	tag "$sampleID"
	label 'bracken'
	label 'long_run'
	publishDir "${params.outdir}/Kraken", mode: 'copy', pattern: "*.bracken",
        saveAs: { filename -> "${sampleID}/${filename}" }


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
	stub:
		bracken_output = sampleID + ".bracken"
		"""
		cat > ${bracken_output} <<EOF
name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads
Bacteria\t2\tD\t1\t0\t1\t1.0
EOF

		cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
		"""
}

process BRACKENMERGE {

	label 'bracken'
	label 'short_run'
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
	stub:
	bracken_merged = "bracken_merged.txt"
	"""
	cat > ${bracken_merged} <<EOF
name\ttaxonomy_id\ttaxonomy_lvl\tstub_sample
Bacteria\t2\tD\t1
EOF

	cat <<-END_VERSIONS > versions.yml
"${task.process}":
  Python: stub
END_VERSIONS
	"""
}