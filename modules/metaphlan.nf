process PREPARE_METAPHLAN {

	executor 'local'
	
	script:

	"""
		cd ${params.metaphlan_db}
		wget https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAAyoJpOgcjop41VIHAGWIVLa/mpa_latest?dl=1
		mv mpa_latest?dl=1 mpa_latest
	
	"""
}

process METAPHLAN {

   label 'metaphlan'
   scratch true
   publishDir "${params.outdir}/${sampleID}/Metaphlan3", mode: 'copy'

   input:
   tuple val(sampleID), file(left_reads), file(right_reads), file(unpaired)

   output:
   path(metaphlan_out), emit: outputMetaphlan
   tuple val(sampleID), file(sam_out), emit: outputMetaphlanBowtie
   tuple val(sampleID), file('v_metaphlan.txt'), emit: version_metaphlan

   script:

   metaphlan_out = sampleID + ".out"
   bowtie_out = sampleID + "bowtie2.txt"
   sam_out = sampleID + ".sam"
   """
     metaphlan --version &> v_metaphlan.txt
     zcat $left_reads > left.fq
     zcat $right_reads > right.fq
     metaphlan left.fq,right.fq --bowtie2db ${params.metaphlan_db} -x mpa_v30_CHOCOPhlAn_201901 --samout $sam_out --bowtie2out $bowtie_out --nproc ${task.cpus} -o $metaphlan_out --input_type fastq
     rm *.fq
   """
}

process ABUNDANCEMERGE {

	publishDir "${params.outdir}/Metaphlan3", mode: 'copy'

	input:
	path(results)

	output:
	file(abundances)

	script:
	abundances = "metaphlan_abundances.txt"

	"""
		merge_metaphlan_tables.py ${results.join(" ")} > $abundances
	"""
}