process FASTQC {

	label 'fastqc'

	publishDir "${params.outdir}/${sampleID}/FastQC", mode: 'copy'

	input:
	tuple val(sampleID), file(left), file(right)

	output:
	path('*_fastqc.{zip,html}')

	script:
	"""
	fastqc --quiet --threads $task.cpus $left $right
	"""
}

process TRIMREADS {

	label 'bbmap'

	scratch true
	input:
	tuple val(sampleID),file(left),file(right)
		
	output:
	tuple val(sampleID),file(left_trimmed),file(right_trimmed), emit: filterPEReads
	tuple val(sampleID),file(unpaired), emit: filterSEReads
	//path bbduk_adapter_stats
			
	script:
	bbduk_adapter_stats = sampleID + ".bbduk.adapter.stats"

	left_trimmed = left.getBaseName() + "_trimmed.fastq.gz"
	right_trimmed = right.getBaseName() + "_trimmed.fastq.gz"

	unpaired = sampleID + "_unpaired.fastq.gz"

	"""
	bbduk.sh stats=$bbduk_adapter_stats threads=${task.cpus} in=${left} in2=${right} out1=${left_trimmed} out2=${right_trimmed} outs=$unpaired ref=${params.adapters} ktrim=r k=23 mink=11 hdist=1 minlength=${params.min_read_length} tpe tbo
	"""
}

process CLEANPEREADS {

	label 'bbmap'

	scratch true

	input:
	tuple val(sampleID),file(left),file(right)

	output:
	tuple val(sampleID),file(left_clean),file(right_clean)

	script:
	left_clean = left.getBaseName() + "_clean.fastq.gz"
	right_clean = right.getBaseName() + "_clean.fastq.gz"
	artifact_stats = sampleID + ".bbduk.artifacts.stats"
		
	"""
	bbduk.sh stats=$artifact_stats threads=${task.cpus} in=${left} in2=${right} k=31 ref=artifacts,phix ordered cardinality out1=${left_clean} out2=${right_clean} minlength=${params.min_read_length}
	"""
}

process CLEANSEREADS {

	label 'bbmap'

	scratch true

	input:
	tuple val(sampleID),file(unpaired)

	output:
	tuple val(sampleID),file(unpaired_clean)

	script:

	unpaired_clean = unpaired.getBaseName() + "_clean.fastq.gz"

	"""
	bbduk.sh threads=${task.cpus} in=${unpaired}  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_clean} minlength=${params.min_read_length}

	"""
}

process FILTERREADS {

	label 'bowtie2'
	scratch true
	//publishDir "${params.outdir}/${sampleID}/reads_clean", mode: 'copy'

	input:
	tuple val(sampleID),file(left),file(right),file(unpaired)
	path(bwt_files)
	each genome
	output:
	tuple val(sampleID),path(left_clean),path(right_clean),path(unpaired_clean), emit: cleaned_reads
	path(bowtie_log), emit: filterlog

	script:
	left_clean = sampleID + ".clean.R1.fastq.gz"
	right_clean = sampleID + ".clean.R2.fastq.gz"
	unpaired_clean = sampleID + ".clean.unpaired.fastq.gz"
	bowtie_log = sampleID + ".txt"
	"""
	bowtie2 -x $genome -1 $left -2 $right -U $unpaired -S /dev/null --no-unal -p ${task.cpus} --un-gz $unpaired_clean  --un-conc-gz ${sampleID}.clean.R%.fastq.gz 2> $bowtie_log
	"""
}

process MULTIQC1 {

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	label 'multiqc'

	input:
	path ('*')
	output:
	file("multiqc_report.html")

	script:

	"""
		cp $params.logo .
                cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .
	"""
}

process MULTIQC2 {

	publishDir "${params.outdir}/MultiQC", mode: 'copy'

	label 'multiqc'

	input:
	path ('*')
	path('*')
	output:
	file("multiqc_report.html")

	script:

	"""
		cp $params.logo .
                cp $baseDir/assets/multiqc_config.yaml multiqc_config.yaml

		multiqc .
	"""
}