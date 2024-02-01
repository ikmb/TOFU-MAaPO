process SYLPH_SKETCH {
scratch params.scratch
tag "$sampleID"
label 'sylph_sketch'

input:
	tuple val(meta), path(reads)

output:
    tuple val(meta), path(sylsp_output), emit: sylph_sketches
	path('versions.yml'), emit: version

script:
	sampleID = meta.id
	sylph_report = sampleID + ".sylph_profile.tbl"
    def args = meta.single_end ? "${reads[0]}.sylsp" : "${reads[0]}.paired.sylsp"
    sylsp_output = sampleID + ".sylsp"
	//sylph_log = sampleID + "_sylph.log"

	if (!meta.single_end) {
	"""

	sylph sketch \
		-1 ${reads[0]} \
		-2 ${reads[1]} \
		-t ${task.cpus} 

    mv $args $sylsp_output

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	sylph: \$(sylph --version 2>&1 | sed -e "s/sylph //g" )
	END_VERSIONS

	"""
	} else {
	"""
	sylph sketch \
		-r ${reads[0]} \
		-t ${task.cpus} 

    mv $args $sylsp_output

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	sylph: \$(sylph --version 2>&1 | sed -e "s/sylph //g" )
	END_VERSIONS
	
	"""	
	}
}

process SYLPH_PROFILING {
scratch params.scratch
tag params.sylph_merge ? "all" : "${meta.id}"
label 'sylph_profile'
publishDir "${params.outdir}/sylph", mode: 'copy', pattern: "*.tbl"

input:
	tuple val(meta), path(sylsp)

output:
    path('*.tbl'), emit: results
	path('versions.yml'), emit: version

script:
    def output_name = params.sylph_merge ? "sylph_merged_profiles.tbl" : "${meta.id}_profile.tbl"
	"""

	sylph profile \
        *.sylsp \
        ${params.sylph_db} \
        -o $output_name \
		-t ${task.cpus} 

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	sylph: \$(sylph --version 2>&1 | sed -e "s/sylph //g" )
	END_VERSIONS

	"""
}