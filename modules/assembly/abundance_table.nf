process abundance_persample{

	label 'bowtie2'
	//scratch params.scratch
    scratch false

	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/MAG_abundance", mode: 'copy'

	input:
        tuple val(meta), file(fafile), file(finalcontigs), path(reads)

	output:
		tuple val(meta), file(fcontigs), file(depthout), emit: maps
		tuple val(meta), file(mappingbam), emit: bam

	script:
        sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"


        refined_bins_file = sampleID + "_refined_bins.fna"

		depthout = sampleID + '_contigs_depth.txt'
		mappingbam = sampleID + '_mapping_final.bam'

		if (!params.single_end) {  
    		"""
            #merge refined bins to one file
            cat refined_bins/*.fa > $refined_bins_file

            mkdir -p bin_coverage

            #create bowtie index from refined_bins
            bowtie2-build $refined_bins_file bin_coverage/bwi_rfbins

            #map qced reads onto contigs
            bowtie2 -p ${task.cpus} \
                    --sensitive \
                    -x bin_coverage/bwi_rfbins \
                    -1 $left_clean \
                    -2 $right_clean \
                    -U $single_clean \
                    -S bin_coverage/${sampleID}.sam

            samtools sort -@ ${task.cpus} \
                          -o bin_coverage/${sampleID}.bam \
                          bin_coverage/${sampleID}.sam

            samtools index bin_coverage/${sampleID}.bam

            rm bin_coverage/${sampleID}.sam

			jgi_summarize_bam_contig_depths bin_coverage/${sampleID}.bam --outputDepth $depthout
			"""
		} else {
			"""
			#merge refined bins to one file
            cat refined_bins/*.fa > $refined_bins_file

            mkdir -p bin_coverage

            #create bowtie index from refined_bins
            bowtie2-build $refined_bins_file bin_coverage/bwi_rfbins

            #map qced reads onto contigs
            bowtie2 -p ${task.cpus} \
                    --sensitive \
                    -x bin_coverage/bwi_rfbins \
                    -U $single_clean \
                    -S bin_coverage/${sampleID}.sam
            
            samtools sort -@ ${task.cpus} \
                          -o bin_coverage/${sampleID}.bam \
                          -n \
                          bin_coverage/${sampleID}.sam            
            
            samtools index bin_coverage/${sampleID}.bam


            rm bin_coverage/${sampleID}.sam
            
			jgi_summarize_bam_contig_depths bin_coverage/${sampleID}.bam --outputDepth $depthout
            """		
		}
}
//rm bin_coverage/${sampleID}.sam
