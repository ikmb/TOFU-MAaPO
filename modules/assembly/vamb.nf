	process VAMB_CATALOGUE {
        label 'vamb'
		scratch params.scratch

        input:
            path(contigs)

        output:
			path(catalogue), emit: catalogue

        script:
			catalogue = "collected_catalogue.fna.gz"

        	"""
        	concatenate.py $catalogue ${contigs.join(" ")} --keepnames
        	"""
    }

	process VAMB_CATALOGUE_INDEX {
        label 'default_highmemory'
		scratch params.scratch

        input:
            path(catalogue)

        output:
			tuple path(catalogue), path(catalogue_index), emit: catalogue

        shell:
			catalogue_index = "catalogue.mmi"

        	"""
			minimap2 -d !{catalogue_index} !{catalogue} -m 2000 --print-qname
			# make index #uses default 3 threads
        	"""
    }
//-I100G
process VAMB_MAPPING{

	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(meta), file(fcontigs), path(reads)
		tuple path(catalogue), path(catalogue_index)

	output:
		path(depthout), emit: counttable 

		val(meta), emit: sampleid
		tuple val(meta), file(fcontigs), file(depthout), emit: maps
		tuple val(meta), file(mappingbam), file(mappingbam_index), emit: bam

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_minimap.bam'
		mappingbam_index = sampleID + '_mapping_minimap.bam.bai'

		sample_total_reads = sampleID + '_totalreads.txt'
		if (!params.single_end) {  
    		"""
				#minimap2 -I100G -d catalogue.mmi $catalogue; # make index
				minimap2 -t ${task.cpus} -N 50 -ax sr  $catalogue_index $left_clean $right_clean | samtools view -F 3584 -b --threads ${task.cpus} | samtools sort > $mappingbam 2> error.log # -n 
				samtools index $mappingbam
				jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			"""
		} else {
			"""	
				#minimap2 -d catalogue.mmi $catalogue; # make index
				minimap2 -t ${task.cpus} -N 5 -ax sr $catalogue_index $single_clean | samtools view -F 3584 -b --threads ${task.cpus}| samtools sort  > $mappingbam 2> error.log #-n
				samtools index $mappingbam
				jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout
			"""		
		}
}

process VAMB_COLLECT_DEPTHS {
	label 'default'
    scratch params.scratch
	//publishDir "${params.outdir}/${sampleID}/vamb", mode: 'copy'

	input:
		path(depthout)

	output:
		path(alldepths), emit: alldepths
	script:

		alldepths = 'all_depths.tsv'

		"""
		Rscript ${baseDir}/bin/collectmapping.R $alldepths
		sed -i "s/[.]var/-var/g" $alldepths
		"""
}

process VAMB {

	label 'vamb'
	scratch params.scratch

	input:
        tuple path(catalogue), path(catalogue_index)
		path(alldepths)

	output:
		path(cluster_table), emit: all_samples_clustertable

	script:
		cluster_table = 'all_vamb_contigs_to_bin.tsv'

		"""
    	vamb --outdir bin --fasta $catalogue --jgi $alldepths -o _${params.contig_sep}_ 
		mv bin/clusters.tsv $cluster_table
		"""
}

//#_k141_ #--bamfiles mappingbam -o C --minfasta 200000

process VAMB_CONTIGS_SELECTION{
	
	label 'default'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/${sampleID}/vamb", mode: 'copy', enabled: params.publish_rawbins
	
	input:
		path(all_cluster_table)
		val(meta)
	output:
		tuple val(meta), file(persample_clustertable), emit: persample_clustertable

	script:
		sampleID = meta.id

		persample_clustertable = sampleID + '_vamb_contigs_to_bin.tsv'

		"""
    	grep $sampleID $all_cluster_table > $persample_clustertable
		"""
}