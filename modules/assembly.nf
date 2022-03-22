
	process MEGAHIT {

	publishDir "${params.outdir}/${sampleID}/Megahit", mode: 'copy'
	scratch params.scratch
	label 'megahit'
	tag "$sampleID"

	input:
		tuple val(meta), path(reads)

	output:
		path('**/*'), emit: outputfolder
		tuple val(sampleID), file("ut-repfix/final.contigs.fa"), path('*_clean.fastq.gz', includeInputs: true), emit: contigs

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

		fq_left = sampleID + "_1.fq"
    	fq_right = sampleID + "_2.fq"
    	fq_single = sampleID + "_single.fq"

		if (!params.single_end) {  
		"""
		zcat $unpaired_clean > $fq_single
		zcat $left_clean > $fq_left
		zcat $right_clean > $fq_right
		megahit -1 $fq_left \
				-2 $fq_right \
				-r $fq_single \
				-m 0.95 \
				-o $sampleID \
				-out-repfix $sampleID \
				-t ${task.cpus}	
		"""
		} else {
		"""	
		zcat $unpaired_clean > $fq_single

		megahit	-r $fq_single \
				-m 0.95 \
				-o $sampleID \
				-out-repfix $sampleID \
				-t ${task.cpus}	
		"""			
		}
	}

	process filtercontigs {
	scratch params.scratch
	tag "$sampleID"

	input:
	tuple val(sampleID), file("ut-repfix/final.contigs.fa"), path(reads)

	output:
	tuple val(sampleID), file("fcontigsfiltered.fa"), path(reads), emit: contigs

	script:
	"""
	python3 ${baseDir}/bin/contigfilterbylen.py 1500 ut-repfix/final.contigs.fa > fcontigsfiltered.fa
	"""
	}

process MAPPING{

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(sampleID), file(fcontigs), path(reads)

	output:
		tuple val(sampleID), file(fcontigs), file(depthout), emit: maps
		tuple val(sampleID), file(mappingbam), emit: counttable

	script:
		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_final.bam'

		if (!params.single_end) {  
    		"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -1 $left_clean -2 $right_clean -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ 5 -o $mappingbam

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout
			"""
		} else {
			"""
			#build and index
			bowtie2-build $fcontigs ${sampleID}_mapping --threads ${task.cpus}
			bowtie2 -p ${task.cpus} -x ${sampleID}_mapping -U $single_clean -S ${sampleID}_mapped.sam |& tee -a ${sampleID}.txt
			samtools view -u ${sampleID}_mapped.sam | samtools sort -m 7G -@ 5 -o $mappingbam

			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout
			"""		
		}
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
	//scratch params.scratch
	scratch false
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