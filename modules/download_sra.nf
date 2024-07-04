process download_sra {
	label 'local_download'
	tag "$sampleID"
	publishDir "${params.outdir}/raw_reads", mode: 'copy', enabled: params.publish_rawreads

	input:
	tuple val(meta), val(reads)

	output:
	tuple val(meta), path("*_raw.fastq.gz") , emit: reads

    beforeScript = """
        sleepFunc() {
            sleep ${params.retryDelay} * delay
        }
        """

	script:
	sampleID = meta.id
    def delay = task.attempt
	leftnewname = sampleID + "_1_raw.fastq.gz"
    rightnewname = sampleID + "_2_raw.fastq.gz"
	unpairednewname = sampleID + "_unpaired_raw.fastq.gz"
    if (meta.single_end) {
		"""
        if [[ "${reads}" = "ftp://ftp."*  ]]; then
            curl -s ${reads} -o $unpairednewname --retry 10 -C -
        else
            curl -s ftp://ftp.sra.ebi.ac.uk${reads} -o $unpairednewname --retry 10 -C -
        fi
        
        """
    } else {
		"""
        if [[ "${reads[0]}" = "ftp://ftp."*  ]]
        then
            curl -s ${reads[0]} -o $leftnewname --retry 10 -C -
            curl -s ${reads[1]} -o $rightnewname --retry 10 -C -
            if [[ -n "${reads[2]}" && "${reads[2]}" != "null" ]]; then 
                curl -s ${reads[2]} -o $unpairednewname --retry 10 -C -
            fi
        else
            curl -s ftp://ftp.sra.ebi.ac.uk${reads[0]} -o $leftnewname --retry 10 -C -
            curl -s ftp://ftp.sra.ebi.ac.uk${reads[1]} -o $rightnewname --retry 10 -C -
            if [[ -n "${reads[2]}" && "${reads[2]}" != "null" ]]; then  
                curl -s ftp://ftp.sra.ebi.ac.uk${reads[2]} -o $unpairednewname --retry 10 -C -
            fi
        fi
        
		"""
    }
    stub:
        sampleID = meta.id
        def delay = task.attempt
        leftnewname = sampleID + "_1_raw.fastq.gz"
        rightnewname = sampleID + "_2_raw.fastq.gz"
        unpairednewname = sampleID + "_unpaired_raw.fastq.gz"
		if (meta.single_end) {
			"""
			touch ${unpairednewname}
			"""
		} else {
			"""
			touch ${leftnewname}
            touch ${rightnewname}
            touch ${unpairednewname}
			"""
		}

}
