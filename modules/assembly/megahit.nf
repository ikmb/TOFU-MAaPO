process MEGAHIT {

	publishDir "${params.outdir}/${sampleID}/Megahit", mode: 'copy'
	scratch params.scratch
	label 'megahit'
	tag "$sampleID"

	input:
		tuple val(meta), path(reads)

	output:
		path('output/*'), emit: outputfolder
		tuple val(meta), file(output_final_contigs), path('*_clean.fastq.gz', includeInputs: true), emit: contigs

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

		fq_left = sampleID + "_1.fq"
    	fq_right = sampleID + "_2.fq"
    	fq_single = sampleID + "_single.fq"

		replacer = ">" + sampleID + "_k"
		output_final_contigs = sampleID + "_final.contigs.fa"

		if (!params.single_end) {  
		    """
		    zcat $unpaired_clean > $fq_single
		    zcat $left_clean > $fq_left
		    zcat $right_clean > $fq_right
            
		    megahit -1 $fq_left \
				-2 $fq_right \
				-r $fq_single \
				-o output \
				-m 0.95 \
				-t ${task.cpus}	

		    rm $fq_single
		    rm $fq_left
		    rm $fq_right

		    gawk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa > $output_final_contigs
		    """
		} else {
		    """	
		    zcat $unpaired_clean > $fq_single

		    megahit	-r $fq_single \
				-m 0.95 \
				-o output \
				-t ${task.cpus}	

		    rm $fq_single

		    gawk '{gsub (/^>k/, "/$replacer");print}' output/final.contigs.fa > $output_final_contigs
		    """			
		}
}