process FAIRY_SKETCH {
	label 'fairy'
	label 'short_run'
	scratch params.scratch
	tag "${meta.id}"

	input:
		tuple val(meta), path(reads)

	output:
		tuple val(meta), path("*.bcsp"), emit: sketch
		path("versions.yml"),            emit: versions

	script:
		def sampleID = meta.id
		//like strobealign for VAMB, only the read pairs are used for paired-end samples
		def read_args = meta.single_end ?
			"-r ${sampleID}_single_clean.fastq.gz" :
			"-1 ${sampleID}_R1_clean.fastq.gz -2 ${sampleID}_R2_clean.fastq.gz"
		"""
		fairy sketch $read_args -S $sampleID -d . -t ${task.cpus}

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		fairy: \$(fairy --version | sed -e "s/fairy //g")
		END_VERSIONS
		"""
	stub:
		def sampleID = meta.id
		"""
		touch ${sampleID}.paired.bcsp
		echo "FAIRY_SKETCH_stub" > versions.yml
		"""
}

process FAIRY_COVERAGE {
	label 'fairy'
	label 'long_run'
	scratch params.scratch
	tag "$key"

	input:
		tuple val(key), path(contigs), path(sketches)

	output:
		tuple val(key), path(depthout), emit: depth    // MetaBAT2 format (jgi_summarize_bam_contig_depths layout)
		path("versions.yml"),           emit: versions

	script:
		depthout = key + '_fairy_depth.txt'
		//sorted sketches give a stable sample column order
		def sketch_files = sketches.collect{ it.toString() }.sort().join(' ')
		"""
		fairy coverage $contigs $sketch_files -t ${task.cpus} -o $depthout

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		fairy: \$(fairy --version | sed -e "s/fairy //g")
		END_VERSIONS
		"""
	stub:
		depthout = key + '_fairy_depth.txt'
		"""
		touch $depthout
		echo "FAIRY_COVERAGE_stub" > versions.yml
		"""
}
