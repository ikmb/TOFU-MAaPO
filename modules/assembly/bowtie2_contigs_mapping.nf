process CONTIGS_MAPPING{

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(sampleID), file(fcontigs), path(reads)

	output:
		tuple val(sampleID), file(fcontigs), file(depthout), emit: maps
		tuple val(sampleID), file(mappingbam), emit: bam

	script:
		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_final.bam'

		if (!params.single_end) {  
    		"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -1 $left_clean -2 $right_clean -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ 5 -o $mappingbam

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout
			"""
		} else {
			"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ 5 -o $mappingbam

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout
			"""		
		}
}