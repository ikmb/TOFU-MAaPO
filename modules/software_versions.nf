process SOFTWARE_VERSIONS {
	label 'default'
	scratch params.scratch
	publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: "*.txt"
	input:
		path('collated_versions.yml')

	output:
		path('TOFU-MAaPO_software_versions.txt')

	script:
		"""
		cat collated_versions.yml | sed -e 's/^[ \t]*//' | sed -n '/END_VERSIONS/!p' > TOFU-MAaPO_software_versions.txt
		"""
}