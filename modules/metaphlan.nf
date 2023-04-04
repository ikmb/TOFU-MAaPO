process PREPARE_METAPHLAN {

	executor 'local'
  label 'local_run'
  output: 
    val 'true', emit: readystate
	script:

	"""
		cd ${params.metaphlan4_db}

    metaphlan --install --force_download --bowtie2db ${params.metaphlan4_db} --nproc ${task.cpus}
    wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest
		##Dropbox mirror
    #wget https://www.dropbox.com/sh/7qze7m7g9fe2xjg/AAAyoJpOgcjop41VIHAGWIVLa/mpa_latest?dl=1
    #mv mpa_latest?dl=1 mpa_latest

    #Zenodo mirror
    #wget https://zenodo.org/record/3957592/files/mpa_latest?download=1
		#mv mpa_latest?download=1 mpa_latest
	
	"""
}

process METAPHLAN {

   label 'metaphlan'
   tag "$sampleID"
   scratch params.scratch
   //TODO: remove bam from output
   //publishDir "${params.outdir}/${sampleID}/Metaphlan4", mode: 'copy'

   input:
    tuple val(meta), path(reads)
    each ready

   output:
    path(metaphlan_out), emit: outputMetaphlan
    tuple val(sampleID), file(bam_out), emit: outputMetaphlanBowtie
    tuple val(sampleID), file('v_metaphlan.txt'), emit: version_metaphlan
   script:
    sampleID = meta.id
    metaphlan_out = sampleID + "_metaphlan.out"
    bowtie_out = sampleID + "_metaphlan_bowtie2.txt"
    sam_out = sampleID + "_metaphlan.sam"
    bam_out = sampleID + "_metaphlan.bam"

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

    phlan_left = sampleID + "_1.fq"
    phlan_right = sampleID + "_2.fq"
    phlan_single = sampleID + "_single.fq"

    if (!params.single_end) {  
      """
      METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
      DEFAULT_DB_FOLDER=${params.metaphlan_db}

      metaphlan --version &> v_metaphlan.txt
     
      zcat ${left_clean} > $phlan_left
      zcat ${right_clean} > $phlan_right
    
      #check if unpaired/single reads are present. 
      # At least 4 reads in "_single_clean.fastq.gz"
      zcat ${unpaired_clean} | head -n 4 &> /dev/null
      if [ $? -eq 0 ];then
        zcat ${unpaired_clean} > $phlan_single

        metaphlan $phlan_left,$phlan_right,$phlan_single \
          --bowtie2db ${params.metaphlan_db} \
          --samout $sam_out \
          --bowtie2out $bowtie_out \
          --stat_q 0.2 \
          --force \
          -t rel_ab_w_read_stats \
          --nproc ${task.cpus} \
          -o $metaphlan_out \
          --input_type fastq \
          --offline

      else
        metaphlan $phlan_left,$phlan_right \
          --bowtie2db ${params.metaphlan_db} \
          --samout $sam_out \
          --bowtie2out $bowtie_out \
          --stat_q 0.2 \
          --force \
          -t rel_ab_w_read_stats \
          --nproc ${task.cpus} \
          -o $metaphlan_out \
          --input_type fastq \
          --offline
      fi
      
      rm *.fq
      samtools view -S -b $sam_out > $bam_out
      rm $sam_out

      """
    } else {
      """
      METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
      DEFAULT_DB_FOLDER=${params.metaphlan_db}

      metaphlan --version &> v_metaphlan.txt
      
      zcat ${unpaired_clean} > $phlan_single

      metaphlan $phlan_single \
        --bowtie2db ${params.metaphlan_db} \
          --samout $sam_out \
          --bowtie2out $bowtie_out \
          --stat_q 0.2 \
          --force \
          -t rel_ab_w_read_stats \
          --nproc ${task.cpus} \
          -o $metaphlan_out \
          --input_type fastq \
          --offline

      rm *.fq
      samtools view -S -b $sam_out > $bam_out
      rm $sam_out

      """
    }
}

process ABUNDANCE_REL_MERGE {
  label 'default'
	publishDir "${params.outdir}/Metaphlan4", mode: 'copy'
  scratch params.scratch

	input:
	  path(results)

	output:
	  file(abundances)

	script:
	  abundances = "metaphlan_rel_abundances.txt"

	  """
    /opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/merge_rel_reads.py ${results.join(" ")} > $abundances
		#merge_metaphlan_tables.py ${results.join(" ")} > $abundances
	  """
}

process ABUNDANCE_ABS_MERGE {
  label 'default'
	publishDir "${params.outdir}/Metaphlan4", mode: 'copy'
  scratch params.scratch

	input:
	  path(results)

	output:
	  file(abundances)

	script:
	  abundances = "metaphlan_abs_abundances.txt"

	  """
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/merge_abs_reads.py ${results.join(" ")} > $abundances
	  """
}
