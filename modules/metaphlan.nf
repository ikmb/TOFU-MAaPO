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
   tag "$sampleID"
   //scratch true
   publishDir "${params.outdir}/${sampleID}/Metaphlan3", mode: 'copy'

   input:
   tuple val(sampleID), file(left_reads), file(right_reads), file(unpaired)

   output:
   path(metaphlan_out), emit: outputMetaphlan
   tuple val(sampleID), file(bam_out), emit: outputMetaphlanBowtie
   tuple val(sampleID), file('v_metaphlan.txt'), emit: version_metaphlan
   //path('*'), emit: metaphlanouts
   script:

   metaphlan_out = sampleID + ".out"
   bowtie_out = sampleID + "_bowtie2.txt"
   sam_out = sampleID + ".sam"
   bam_out = sampleID + ".bam"
   """
     metaphlan --version &> v_metaphlan.txt
     zcat $left_reads > left.fq
     zcat $right_reads > right.fq
     metaphlan left.fq,right.fq --bowtie2db ${params.metaphlan_db} -x mpa_v30_CHOCOPhlAn_201901 --samout $sam_out --bowtie2out $bowtie_out --stat_q 0.2 --force -t rel_ab_w_read_stats --nproc ${task.cpus} -o $metaphlan_out --input_type fastq
     rm *.fq
     samtools view -S -b $sam_out > $bam_out
     rm $sam_out
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