process DOWNLOAD_HUMANGENOME {
	executor 'local'
    label 'local_run'
	label 'short_run'
    output: 
        val 'true', emit: readystate
	script:

	"""
    if [ ! -d ${params.genomes.human.bowtie_index} ]; then
		mkdir -p ${params.genomes.human.bowtie_index};
    fi
	cd ${params.genomes.human.bowtie_index}

    wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_decoy_as.zip
	gunzip GRCh38_noalt_decoy_as.zip
	"""
}