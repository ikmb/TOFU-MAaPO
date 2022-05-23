
	process MEGAHIT {

	publishDir "${params.outdir}/${sampleID}/Megahit", mode: 'copy'
	scratch params.scratch
	label 'megahit'
	tag "$sampleID"

	input:
		tuple val(meta), path(reads)

	output:
		path('output/*'), emit: outputfolder
		tuple val(sampleID), file(output_final_contigs), path('*_clean.fastq.gz', includeInputs: true), emit: contigs

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		unpaired_clean = sampleID + "_single_clean.fastq.gz"

		fq_left = sampleID + "_1.fq"
    	fq_right = sampleID + "_2.fq"
    	fq_single = sampleID + "_single.fq"

		replacer = ">" + sampleID + "_k"
		output_final_contigs = sampleID + "_final.contigs.fa"
		if (!params.single_end) {  
		"""
		zcat $unpaired_clean > $fq_single
		zcat $left_clean > $fq_left
		zcat $right_clean > $fq_right
		megahit -1 $fq_left \
				-2 $fq_right \
				-r $fq_single \
				-o output \
				-m 0.95 \
				-t ${task.cpus}	

		rm $fq_single
		rm $fq_left
		rm $fq_right

		gawk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa > $output_final_contigs
		"""
		} else {
		"""	
		zcat $unpaired_clean > $fq_single

		megahit	-r $fq_single \
				-m 0.95 \
				-o output \
				-t ${task.cpus}	

		rm $fq_single

		gawk '{gsub (/^>k/, "/$replacer");print}' output/final.contigs.fa > $output_final_contigs
		"""			
		}
	}

	process filtercontigs {
	scratch params.scratch
	tag "$sampleID"

	input:
	tuple val(sampleID), file(finalcontigs), path(reads)

	output:
	tuple val(sampleID), file("fcontigsfiltered.fa"), path(reads), emit: contigs

	script:
	"""
	python3 ${baseDir}/bin/contigfilterbylen.py 1500 $finalcontigs > fcontigsfiltered.fa
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
	export GTDBTK_DATA_PATH=${params.gtdbtk_reference}
	gtdbtk classify_wf --cpus ${task.cpus} --genome_dir . --extension fa --out_dir all.bins.gtdbtk_output --pplacer_cpus 1

	gawk -F "\t" '{ sub(/.*;s__/, "s__", \$2); print \$1 "\t" \$2 }' all.bins.gtdbtk_output/gtdbtk.bac120.summary.tsv > all.bins.gtdbtk_output/parsed_bac120_summary.tsv
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