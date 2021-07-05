

	process MEGAHIT {

	publishDir "${OUTDIR}/${sampleID}/Megahit", mode: 'copy'
	//scratch true
	label 'megahit'

	input:
	tuple val(sampleID),path(left),path(right),path(unpaired)
	output:
	file("**/*"), emit: outputfolder
	tuple val(sampleID), file("ut-repfix/final.contigs.fa"), path(left),path(right),path(unpaired), emit: contigs
	script:
	"""
	zcat $unpaired > unpaired.fq
	zcat $left > left.fq
    zcat $right > right.fq
	megahit -1 left.fq -2 right.fq -r unpaired.fq -m 0.95 -out-repfix $sampleID -t ${task.cpus}	
	"""
	}

	process MAPPING {

	label 'bowtie2'
	scratch true
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
	scratch true
	publishDir "${OUTDIR}/${sampleID}/Metabat2", mode: 'copy'

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

	publishDir "${OUTDIR}/${sampleID}/Metabat2", mode: 'copy'
	scratch true
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
	//scratch true
	publishDir "${OUTDIR}/${sampleID}/checkm", mode: 'copy'

	input: 
	tuple val(sampleID), file(fafile)

	output:
	
	shell:
	"""
	checkm lineage_wf -t ${task.cpus} -x fa . ./bins
	"""
	}

	process GTDBTK {

	label 'gtdbtk'
	//scratch true
	publishDir "${OUTDIR}/${sampleID}/gtdbtk", mode: 'copy'

	input: 
	tuple val(sampleID), file(fafile)

	output:
	
	shell:
	"""
	export GTDBTK_DATA_PATH="${baseDir}/databases/release202"
	gtdbtk classify_wf --cpus ${task.cpus} --genome_dir . --extension fa --out_dir all.bins.gtdbtk_output --pplacer_cpus 1
	"""
	}

	process getCountTable {
	publishDir "${OUTDIR}/${sampleID}/counttable", mode: 'copy'

	input:
	tuple val(sampleID), file(finalbam)
	output:

	shell:
	"""
	samtools idxstats $finalbam > ${sampleID}_idxstats.txt
    python ${baseDir}/bin/get_count_table.py ${sampleID}_idxstats.txt > counts_${sampleID}.txt
	"""
	}