process download_sra {
	label 'local_download'
	tag "$sampleID"
	publishDir "${params.outdir}/raw_reads", mode: 'copy', enabled: params.publish_rawreads

	input:
	tuple val(meta), val(reads)

	output:
	tuple val(meta), path("*_raw.fastq.gz") , emit: reads

	script:
	sampleID = meta.id

	leftnewname = sampleID + "_1_raw.fastq.gz"
    rightnewname = sampleID + "_2_raw.fastq.gz"
	unpairednewname = sampleID + "_unpaired_raw.fastq.gz"
    if (meta.single_end) {
		"""
        if [[ "${reads[0]}" = "ftp://ftp."*  ]]; then
            wget --quiet -c ${reads[0]} -O $unpairednewname
        else
            wget --quiet -c ftp://ftp.sra.ebi.ac.uk${reads[0]} -O $unpairednewname
        fi
        
        """
    } else {
		"""
        if [[ "${reads[0]}" = "ftp://ftp."*  ]]
        then
            wget ${reads[0]} -O $leftnewname
            wget ${reads[1]} -O $rightnewname
            if [[ -n "${reads[2]}" && "${reads[2]}" != "null" ]]; then 
                wget ${reads[2]} -O $unpairednewname
            fi
        else
            wget ftp://ftp.sra.ebi.ac.uk${reads[0]} -O $leftnewname
            wget ftp://ftp.sra.ebi.ac.uk${reads[1]} -O $rightnewname
            if [[ -n "${reads[2]}" && "${reads[2]}" != "null" ]]; then  
                wget ftp://ftp.sra.ebi.ac.uk${reads[2]} -O $unpairednewname
            fi
        fi
        
		"""
    }
}