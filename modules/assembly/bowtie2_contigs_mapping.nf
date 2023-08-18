process CONTIGS_MAPPING{

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), path(reads)

	output:
		tuple val(meta), file(fcontigs), file(depthout), emit: maps
		tuple val(meta), file(mappingbam), emit: bam
		tuple val(meta), file(depthout), file(sample_total_reads), emit: sample_depth
		path("versions.yml"),	optional: true, emit: versions

	script:
		sampleID = meta.id
		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_final.bam'
		sample_total_reads = sampleID + '_totalreads.txt'

		if (!meta.single_end) {  
			"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -1 $left_clean -2 $right_clean -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ ${task.cpus} -o $mappingbam

			grep "reads; of these:" ${sampleID}.txt | sed 's/ reads; of these://' > $sample_total_reads

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
      		bowtie2: \$(bowtie2 --version | awk 'FNR==1' |sed 's/.* //')
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}')
			END_VERSIONS
			"""
		} else {
			"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ ${task.cpus} -o $mappingbam

			grep "reads; of these:" ${sampleID}.txt | sed 's/ reads; of these://' > $sample_total_reads

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
      		bowtie2: \$(bowtie2 --version | awk 'FNR==1' |sed 's/.* //')
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}')
			END_VERSIONS
			"""		
		}
}