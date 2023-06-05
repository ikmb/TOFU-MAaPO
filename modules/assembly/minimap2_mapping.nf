	process MINIMAP2_CATALOGUE {
        label 'vamb'
		scratch params.scratch
		tag "$sampleID"
        input:
		    tuple val(meta), file(fcontigs), path(reads)
        output:
			tuple val(meta), path(catalogue), optional: true, emit: catalogue
			path("versions.yml"), emit: versions

        script:
			sampleID = meta.id
			catalogue = "collected_catalogue.fna.gz"

        	"""
        	concatenate.py $catalogue ${fcontigs.join(" ")} --keepnames

			cat <<-END_VERSIONS > versions.yml
    		"${task.process}":
      		Python: \$(python --version | sed -e "s/Python //g" )
    		END_VERSIONS
        	"""
    }

	process MINIMAP2_CATALOGUE_INDEX {
        label 'default'
		scratch params.scratch
		tag "$sampleID"
        input:
            tuple val(meta), path(catalogue)

        output:
			tuple val(meta), path(catalogue), path(catalogue_index), emit: catalogue
			path("versions.yml"), emit: versions

        script:
			sampleID = meta.id
			catalogue_index = "catalogue.mmi"

        	"""
			minimap2 -I100G -d $catalogue_index $catalogue -m 2000 # make index

			cat <<-END_VERSIONS > versions.yml
    		"${task.process}":
      		minimap2: \$(minimap2 --version)
    		END_VERSIONS
        	"""
    }
process MINIMAP2_MAPPING{

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), path(reads), path(catalogue), path(catalogue_index)

	output:
		path(depthout), 											emit: counttable
		val(meta), 													emit: sampleid
        tuple val(meta), file(depthout), 							emit: sample_depth
		tuple val(meta), file(fcontigs), file(depthout), 			emit: maps
		tuple val(meta), file(mappingbam), file(mappingbam_index), 	emit: bam
		path("error.log"),    										emit: errorlog
		path("versions.yml"), emit: versions

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_minimap.bam'
		mappingbam_index = sampleID + '_mapping_minimap.bam.bai'

		sample_total_reads = sampleID + '_totalreads.txt'
		if (!params.single_end) {  
    		"""
			#minimap2 -I100G -d catalogue.mmi $catalogue; # make index
			minimap2 -t ${task.cpus} -N 50 -ax sr  $catalogue_index $left_clean $right_clean | samtools view -F 3584 -b --threads ${task.cpus} | samtools sort > $mappingbam 2> error.log # -n 
			samtools index $mappingbam
			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
    		"${task.process}":
      		minimap2: \$(minimap2 --version)
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}' )
    		END_VERSIONS

			"""
		} else {
			"""	
			#minimap2 -d catalogue.mmi $catalogue; # make index
			minimap2 -t ${task.cpus} -N 5 -ax sr $catalogue_index $single_clean | samtools view -F 3584 -b --threads ${task.cpus}| samtools sort  > $mappingbam 2> error.log #-n
			samtools index $mappingbam
			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
    		"${task.process}":
      		minimap2: \$(minimap2 --version)
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}')
    		END_VERSIONS
			"""		
		}
}