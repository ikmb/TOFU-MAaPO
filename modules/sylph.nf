process SYLPH_SKETCH {
tag "$sampleID"
label 'sylph_sketch'

input:
	tuple val(meta), path(reads)

output:
    tuple val(meta), path('*.sylsp'), optional: true, emit: sylph_sketches
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	sylph_report = sampleID + ".sylph_profile.tbl"
    def args = meta.single_end ? "${reads[0]}" : "${reads[0]}.paired"
    sylsp_output = args + ".sylsp"
    def runspec = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

	"""
	sylph sketch $runspec -t ${task.cpus} 

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	sylph: \$(sylph --version 2>&1 | sed -e "s/sylph //g" )
	END_VERSIONS

	"""

}

process SYLPH_PROFILING {
tag params.sylph_merge ? "all" : "${meta.id}"
label 'sylph_profile'
publishDir "${params.outdir}/sylph", mode: 'copy', pattern: "*.tbl"

input:
	tuple val(meta), path(sylsp), path(sylph_db)

output:
    path('*.tbl'), optional: true, emit: results
	path('versions.yml'), emit: version

script:
    def output_name = params.sylph_merge ? "sylph_merged_profiles.tbl" : "${meta.id}_profile.tbl"
	"""

	sylph profile \
		-u --read-seq-id 99.1 \
        ${sylph_db} \
		${sylsp.join(" ")} \
        -o $output_name \
		-t ${task.cpus} 

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	sylph: \$(sylph --version 2>&1 | sed -e "s/sylph //g" )
	END_VERSIONS

	"""
}