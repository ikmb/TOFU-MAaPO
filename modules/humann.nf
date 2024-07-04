	process PREPARE_HUMANN {

	label 'local_run'
	output: 
		val 'true', emit: readystate
	script:

	"""
	if [ ! -d ${params.humann_db} ]; then
		mkdir -p ${params.humann_db};
	fi
	cd ${params.humann_db}

	humann_databases --update-config no --download uniref uniref90_diamond ${params.humann_db}
	humann_databases --update-config no --download chocophlan full ${params.humann_db}
	"""
}
	

	process HUMANN {

	label 'humann'
	tag "$sampleID"
	scratch params.scratch
	errorStrategy { (task.exitStatus in [143,137,104,134,139,1] && task.attempt <= maxRetries)  ? 'retry' : 'ignore' }

	input:
		tuple val(meta),path(reads)
		each readymetaphlan
		each readyhumann

	output:
		path('*_genefamilies.tsv'),    optional: true, emit: genefamilies
		path('*_pathabundance.tsv'),   optional: true, emit: pathabundance
		path('*_pathcoverage.tsv'),    optional: true, emit: pathcoverage
		path('*'),                     optional: true, emit: humannouts
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

		phlan_left = sampleID + "_1.fq"
		phlan_right = sampleID + "_2.fq"
		phlan_single = sampleID + "_single.fq"

		genefamilies = sampleID + '_genefamilies.tsv'
		pathabundance = sampleID + '_pathabundance.tsv'
		pathcoverage = sampleID + '_pathcoverage.tsv'
		merged = sampleID + '.fq'

		if (!meta.single_end) {
		"""
		zcat ${left_clean} > $phlan_left
		zcat ${right_clean} > $phlan_right
		
		#check if unpaired/single reads are present
		if [ -s ${unpaired_clean} ];then

			zcat ${unpaired_clean} > $phlan_single
			cat $phlan_left $phlan_right $phlan_single > $merged
		else
			cat $phlan_left $phlan_right > $merged
		fi

		METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
		DEFAULT_DB_FOLDER=${params.metaphlan_db}

		humann --input $merged \
			--output . \
			--remove-temp-output \
			--threads ${task.cpus} \
			--nucleotide-database ${params.humann_db}/chocophlan \
			--protein-database ${params.humann_db}/uniref \
			--metaphlan-options "--bowtie2db ${params.metaphlan_db} --offline" \
			--output-basename $sampleID
		rm *.fq

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		humann: \$(humann --version 2>&1 | sed -e "s/humann v//g")
		END_VERSIONS

		"""

		} else {
		"""
		zcat ${unpaired_clean} > $phlan_single
		
		METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
		DEFAULT_DB_FOLDER=${params.metaphlan_db}

		humann --input $phlan_single \
			--output . \
			--output-basename $sampleID \
			--remove-temp-output \
			--threads ${task.cpus} \
			--nucleotide-database ${params.humann_db}/chocophlan \
			--protein-database ${params.humann_db}/uniref \
			--metaphlan-options "--bowtie2db ${params.metaphlan_db} --offline"
		rm *.fq

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		humann: \$(humann --version 2>&1 | sed -e "s/humann v//g")
		END_VERSIONS

		"""
		}
}

process JOINgenefamilies {
	
	label 'humann'
	publishDir "${params.outdir}/humann", mode: 'copy', pattern: "*.tsv"
	scratch params.scratch
	
	input:
	path('*')

	output:
		path(mergedtable),     emit: abundances
		path("versions.yml"),  emit: versions

	script:
		mergedtable = "humann_merged_genefamilies.tsv"

		"""
		humann_join_tables --input . --output $mergedtable

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS

		"""
}

process JOINpathabundance {

	label 'humann'
	publishDir "${params.outdir}/humann", mode: 'copy', pattern: "*.tsv"
	scratch params.scratch
	
	input:
	path('*')

	output:
		path(mergedtable),     emit: abundances
		path("versions.yml"),  emit: versions

	script:
		mergedtable = "humann_merged_pathabundance.tsv"

		"""
		humann_join_tables --input . --output $mergedtable
		
		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS

		"""
}

process JOINpathcoverage {

	label 'humann'
	publishDir "${params.outdir}/humann", mode: 'copy', pattern: "*.tsv"
	scratch params.scratch

	input:
	path('*')

	output:
		path(mergedtable),     emit: abundances
		path("versions.yml"),  emit: versions

	script:
		mergedtable = "humann_merged_pathcoverage.tsv"

		"""
		humann_join_tables --input . --output $mergedtable

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS
		
		"""
}       