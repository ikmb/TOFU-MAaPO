process PREPARE_METAPHLAN {

	executor 'local'
	
	script:

	"""
		cd ${params.metaphlan_db}
		##Dropbox mirror is currently down
    #wget https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAAyoJpOgcjop41VIHAGWIVLa/mpa_latest?dl=1
    #mv mpa_latest?dl=1 mpa_latest

    #Zenodo mirror
    wget https://zenodo.org/record/3957592/files/mpa_latest?download=1
		mv mpa_latest?download=1 mpa_latest
	
	"""
}

process METAPHLAN {

   label 'metaphlan'
   tag "$sampleID"
   //scratch true
   publishDir "${params.outdir}/${sampleID}/Metaphlan3", mode: 'copy'

   input:
   //tuple val(sampleID), file(left_reads), file(right_reads), file(unpaired)
   tuple val(meta), path(reads)

   output:
    path(metaphlan_out), emit: outputMetaphlan
    tuple val(sampleID), file(bam_out), emit: outputMetaphlanBowtie
    tuple val(sampleID), file('v_metaphlan.txt'), emit: version_metaphlan
    //path('*'), emit: metaphlanouts
   script:
    sampleID = meta.id
    metaphlan_out = sampleID + ".out"
    bowtie_out = sampleID + "_bowtie2.txt"
    sam_out = sampleID + ".sam"
    bam_out = sampleID + ".bam"

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

    phlan_left = sampleID + "_1.fq"
    phlan_right = sampleID + "_2.fq"
    phlan_single = sampleID + "_single.fq"

    if (!params.single_end) {  
      """
      metaphlan --version &> v_metaphlan.txt
     
      zcat ${left_clean} > $phlan_left
      zcat ${right_clean} > $phlan_right
    
      #check if unpaired/single reads are present
      if [ -s ${unpaired_clean} ];then
        zcat ${unpaired_clean} > $phlan_single

        metaphlan $phlan_left,$phlan_right,$phlan_single \
          --bowtie2db ${params.metaphlan_db} \
          -x mpa_v30_CHOCOPhlAn_201901 \
          --samout $sam_out \
          --bowtie2out $bowtie_out \
          --stat_q 0.2 \
          --force \
          -t rel_ab_w_read_stats \
          --nproc ${task.cpus} \
          -o $metaphlan_out \
          --input_type fastq
      else
        metaphlan $phlan_left,$phlan_right \
          --bowtie2db ${params.metaphlan_db} \
          -x mpa_v30_CHOCOPhlAn_201901 \
          --samout $sam_out \
          --bowtie2out $bowtie_out \
          --stat_q 0.2 \
          --force \
          -t rel_ab_w_read_stats \
          --nproc ${task.cpus} \
          -o $metaphlan_out \
          --input_type fastq
      fi
      
      rm *.fq
      samtools view -S -b $sam_out > $bam_out
      rm $sam_out
      """
    } else {
      """
      metaphlan --version &> v_metaphlan.txt
      
      zcat ${unpaired_clean} > $phlan_single

      metaphlan $phlan_single \
        --bowtie2db ${params.metaphlan_db} \
        -x mpa_v30_CHOCOPhlAn_201901 \
        --samout $sam_out \
        --bowtie2out $bowtie_out \
        --stat_q 0.2 \
        --force \
        -t rel_ab_w_read_stats \
        --nproc ${task.cpus} \
        -o $metaphlan_out \
        --input_type fastq
        
      rm *.fq
      samtools view -S -b $sam_out > $bam_out
      rm $sam_out
      """
    }
}

process ABUNDANCE_REL_MERGE {

	publishDir "${params.outdir}/Metaphlan3", mode: 'copy'

	input:
	path(results)

	output:
	file(abundances)

	script:
	abundances = "metaphlan_rel_abundances.txt"

	"""
		merge_metaphlan_tables.py ${results.join(" ")} > $abundances
	"""
}

process ABUNDANCE_ABS_MERGE {

	publishDir "${params.outdir}/Metaphlan3", mode: 'copy'

	input:
	path(results)

	output:
	file(abundances)

	script:
	abundances = "metaphlan_abs_abundances.txt"

	"""
		python3 ${baseDir}/bin/merge_abs_reads.py ${results.join(" ")} > $abundances
	"""
}