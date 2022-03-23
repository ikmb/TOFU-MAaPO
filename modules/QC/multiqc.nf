process MULTIQC1 {

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	label 'multiqc'

	input:
		path('*')
		path('*')
	output:
		file("multiqc_report.html")

	script:

		"""
		cp $params.logo .
        cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .
		"""
}

process MULTIQC2 {

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	label 'multiqc'

	input:
		path('*')
		path('*')
		path('*')
	output:
		file("multiqc_report.html")

	script:

		"""
		cp $params.logo .
        cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .
		"""
}