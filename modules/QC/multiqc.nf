process MULTIQC {

	publishDir "${params.outdir}/MultiQC", mode: 'copy', pattern: "*.html"
	scratch params.scratch
	label 'multiqc'

	input:
		path('*')
	output:
		path("multiqc_report.html")
		path("versions.yml"), emit: versions

	script:

		"""
		cp $params.logo .
        cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	multiqc: \$(multiqc --version| sed -e "s/multiqc, version //g")
		END_VERSIONS
		"""
}
