	process MEGAHIT {

	publishDir "${params.outdir}/${sampleID}/Megahit", mode: 'copy'
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
	megahit -1 left.fq -2 right.fq -r unpaired.fq -m 0.95 -o $sampleID -out-repfix $sampleID -t ${task.cpus}	
	"""
	}

	process filtercontigs {
	scratch params.scratch
	tag "$sampleID"

	input:
	tuple val(sampleID), file("ut-repfix/final.contigs.fa"), path(left),path(right),path(unpaired)

	output:
	tuple val(sampleID), file("fcontigsfiltered.fa"), path(left),path(right),path(unpaired), emit: contigs

	script:
	"""
	python3 ${baseDir}/bin/contigfilterbylen.py 1500 ut-repfix/final.contigs.fa > fcontigsfiltered.fa
	"""
	}
	
	process MAPPING {

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"

	input:
	tuple val(sampleID), file(fcontigs), path(left),path(right),path(unpaired)

	output:
	tuple val(sampleID), file(fcontigs), file(depthout), emit: maps
	tuple val(sampleID), file("${sampleID}_final.bam"), emit: counttable

	script:
	depthout = sampleID + '_depth.txt'
    """
	#build and index
	bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
	bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -1 $left -2 $right -U $unpaired -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
	samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ 5 -o ${sampleID}_final.bam  

	jgi_summarize_bam_contig_depths ${sampleID}_final.bam --outputDepth ${sampleID}_depth.txt
	"""

	}

	process METABAT {

	label 'metabat2'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/Metabat2", mode: 'copy'

	input: 
	tuple val(sampleID), file(fcontigs), file(depthout)

	output:
	tuple val(sampleID), file("${sampleID}_bin.*.fa")

	script:
	"""
	metabat2 -i $fcontigs -a ${sampleID}_depth.txt -o ${sampleID}_bin -t ${task.cpus}
	"""
	}

	process contigs_to_bins {

	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/Metabat2", mode: 'copy'
	
	input: 
	tuple val(sampleID), file(fafile)

	output:
	file("${sampleID}.contigs_to_bins.tsv")

	shell:
	"""
	grep '>' !{fafile} | tr '>' '\t' | tr -d ':' > !{sampleID}.contigs_to_bins.tsv
	"""
	}

	process checkm_all_bins {

	label 'checkm'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/checkm", mode: 'copy'

	input: 
	tuple val(sampleID), file(fafile)

	output:
	file('bins/*')

	shell:
	"""
	checkm lineage_wf -t ${task.cpus} -x fa . ./bins
	"""
	}

	process GTDBTK {

	label 'gtdbtk'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/gtdbtk", mode: 'copy'

	input: 
	tuple val(sampleID), file(fafile)

	output:
	file("all.bins.gtdbtk_output/*")
	shell:
	"""
	export GTDBTK_DATA_PATH=${params.GTDBTKreference}
	gtdbtk classify_wf --cpus ${task.cpus} --genome_dir . --extension fa --out_dir all.bins.gtdbtk_output --pplacer_cpus 1
	"""
	}

	process getCountTable {
	
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/counttable", mode: 'copy'

	input:
	tuple val(sampleID), file(finalbam)

	output:
	file("*.txt")

	shell:
	"""
	samtools idxstats $finalbam > ${sampleID}_idxstats.txt
    python ${baseDir}/bin/get_count_table.py ${sampleID}_idxstats.txt > counts_${sampleID}.txt
	"""
	}