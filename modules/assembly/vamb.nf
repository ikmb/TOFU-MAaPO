	process VAMB_CONCATENATE {
        label 'bowtie2'
		//scratch params.scratch
		scratch false
        input:
            path(contigs)

        output:
			tuple path(catalogue), path(catalogue_index), emit: catalogue

        script:
			catalogue = "collected_catalogue.fna.gz"
			catalogue_index = "catalogue.mmi"

        """
        concatenate.py $catalogue ${contigs.join(" ")} --keepnames
		minimap2 -I100G -d $catalogue_index $catalogue -m 2000 # make index
        """
    }

process VAMB_MAPPING{

	label 'bowtie2'
	//scratch params.scratch
	scratch false
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(sampleID), file(fcontigs), path(reads)
		tuple path(catalogue), path(catalogue_index)

	output:
		path(depthout), emit: counttable 

		val(sampleID), emit: sampleid

	script:
		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_minimap.bam'
		mappingbam_index = sampleID + '_mapping_minimap.bam.bai'

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

    scratch false
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
	//scratch params.scratch
    scratch false
	//publishDir "${params.outdir}/${sampleID}/vamb", mode: 'copy'

	input:
        tuple path(catalogue), path(catalogue_index)
		path(alldepths) //, file(mappingbam), file(mappingbam_index)
	    //tuple val(sampleID), path(catalogue), path(mappingbam)//, path(mappingbam_index)
        //tuple val(sampleID), file(fcontigs), file(depthout)

	output:
	    //tuple val(sampleID), file("${sampleID}_bin/bins/*.fna"), emit: fna_output
		path(cluster_table), emit: all_samples_clustertable

	script:
		cluster_table = 'all_vamb_contigs_to_bin.tsv'

	"""
    	vamb --outdir bin --fasta $catalogue --jgi $alldepths -o _k141_ #--bamfiles mappingbam -o C --minfasta 200000
		mv bin/clusters.tsv $cluster_table
	"""
}



process VAMB_CONTIGS_SELECTION{
	scratch params.scratch
	tag "$sampleID"

	input:
		path(all_cluster_table)
		val(sampleID)
	output:
		tuple val(sampleID), file(persample_clustertable), emit: persample_clustertable

	script:
		persample_clustertable = sampleID + '_vamb_contigs_to_bin.tsv'

	"""
    	grep $sampleID $all_cluster_table > $persample_clustertable
	"""
}