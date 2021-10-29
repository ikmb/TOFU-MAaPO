process metaspades {

	publishDir "${params.outdir}/${sampleID}/metaspades", mode: 'copy'
	scratch params.scratch
	label 'megahit'
	tag "$sampleID"

	input:
	tuple val(sampleID),path(left),path(right),path(unpaired)

	output:
	path("**/*"), emit: outputfolder
	tuple val(sampleID), file("ut-repfix/final.contigs.fa"), path(left),path(right),path(unpaired), emit: contigs

	script:
	"""
	zcat $unpaired > unpaired.fq
	zcat $left > left.fq
	zcat $right > right.fq
	spades.py --meta --pe-1 left.fq --pe-2 right.fq --pe-s unpaired.fq -k 21,33,55 -o ${sample}_spades_out -t ${task_cpus}

	"""
	}