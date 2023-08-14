process PREPARE_METAPHLAN {

	executor 'local'
  label 'local_run'
  output: 
    val 'true', emit: readystate
	script:

	"""
    if [ ! -d ${params.metaphlan_db} ]; then
      mkdir -p ${params.metaphlan_db};
    fi
    
		cd ${params.metaphlan_db}

    #within the humann config the latest version of the metaphlanDB to work is hardcoded to "mpa_vJan21_CHOCOPhlAnSGB_202103", so we need to fake it to be the latest available version
    metaphlan --install --force_download --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db ${params.metaphlan_db} --nproc ${task.cpus}
    echo "mpa_vJan21_CHOCOPhlAnSGB_202103" > mpa_latest

    #wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_latest
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
    path(metaphlan_out),                          optional: true, emit: outputMetaphlan
    tuple val(sampleID), file("*_metaphlan.*am"), optional: true, emit: outputMetaphlanBowtie
    path("versions.yml"),                         optional: true, emit: versions
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

    if (!meta.single_end) {  
      """
      METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
      DEFAULT_DB_FOLDER=${params.metaphlan_db}
    
      zcat ${left_clean} > $phlan_left
      zcat ${right_clean} > $phlan_right
    
      #check if unpaired/single reads are present
      if [[ -f ${unpaired_clean} && \$(zcat ${unpaired_clean} | wc -l) -ge 4 ]]; then
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

      set +e
      {
      samtools view -S -b $sam_out > $bam_out
      if [[ -f ${bam_out} && \$(zcat ${bam_out} | wc -l) -ge 4 ]]; then
        rm $sam_out
      else
        rm $bam_out
      fi
      }
      set -e

      if [ ! -f ${bam_out} ]; then
        >&2 echo "Warning: ${bam_out} could not be created, please check the file ${sam_out}!"
      fi

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
      metaphlan4: \$(metaphlan --version 2>&1 | awk '{print \$3}')
      END_VERSIONS

      """
    } else {
      """
      METAPHLAN_BOWTIE2_DB=${params.metaphlan_db}
      DEFAULT_DB_FOLDER=${params.metaphlan_db}
      
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

      set +e
      {
      samtools view -S -b $sam_out > $bam_out
      if [[ -f ${bam_out} && \$(zcat ${bam_out} | wc -l) -ge 4 ]]; then
        rm $sam_out
      else
        rm $bam_out
      fi
      }
      set -e

      if [ ! -f ${bam_out} ]; then
        >&2 echo "Warning: ${bam_out} could not be created, please check the file ${sam_out}!"
      fi

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
      metaphlan4: \$(metaphlan --version 2>&1 | awk '{print \$3}')
      END_VERSIONS

      """
    }
}

process ABUNDANCE_REL_MERGE {
  label 'default'
	publishDir "${params.outdir}/Metaphlan4", mode: 'copy', pattern: "*.txt"
  scratch params.scratch

	input:
	  path(results)

	output:
	  path(abundances),     emit: abundances
    path("versions.yml"), emit: versions

	script:
	  abundances = "metaphlan_rel_abundances.txt"

	  """
    /opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/merge_rel_reads.py ${results.join(" ")} > $abundances
		#merge_metaphlan_tables.py ${results.join(" ")} > $abundances
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    Python: \$(python --version | sed -e "s/Python //g" )
    END_VERSIONS

	  """
}

process ABUNDANCE_ABS_MERGE {
  label 'default'
	publishDir "${params.outdir}/Metaphlan4", mode: 'copy', pattern: "*.txt"
  scratch params.scratch

	input:
	  path(results)

	output:
	  path(abundances),     emit: abundances
    path("versions.yml"), emit: versions

	script:
	  abundances = "metaphlan_abs_abundances.txt"

	  """
		/opt/conda/envs/ikmb-metagenome-1.2/bin/python3 ${baseDir}/bin/merge_abs_reads.py ${results.join(" ")} > $abundances

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    Python: \$(python --version | sed -e "s/Python //g" )
    END_VERSIONS
    
	  """
}
