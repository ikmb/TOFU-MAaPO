	process VAMB_CATALOGUE {
        label 'vamb'
		scratch params.scratch
		tag "$vamb_key"

        input:
            tuple val(vamb_key), val(contigs)

        output:
			tuple val(vamb_key), path(catalogue), emit: catalogue
			path("versions.yml"), emit: versions
        script:
			catalogue = vamb_key + "_collected_catalogue.fna.gz"

			"""
			concatenate.py $catalogue ${contigs} --keepnames

			cat <<-END_VERSIONS> versions.yml
			"${task.process}":
			Python: \$(python --version | sed -e "s/Python //g" )
			END_VERSIONS

			"""
    }

	process VAMB_CATALOGUE_INDEX {
        label 'default_highmemory'
		cache 'lenient'
		scratch params.scratch
		tag "$vamb_key"

        input:
            tuple val(vamb_key), path(catalogue)

        output:
			tuple  path(catalogue), val(vamb_key), path(catalogue_index), emit: catalogue
			tuple  val(vamb_key), path(catalogue), path(catalogue_index), emit: catalogue_indexfirst
			path("versions.yml"), emit: versions
        shell:
			catalogue_index = "catalogue.mmi"

			"""
			minimap2 -d !{catalogue_index} !{catalogue} -m 2000 --print-qname -I!{params.minimap_indexsize}g
			# make index #uses default 3 threads

			cat <<-END_VERSIONS> versions.yml
			"${task.process}":
			minimap2: \$(minimap2 --version)
			END_VERSIONS

			"""
    }
//-I100G
process VAMB_MAPPING{
	cache 'lenient'
	label 'bowtie2'
	scratch params.scratch
	tag "$sampleID"
	//publishDir "${params.outdir}/${sampleID}/Mapping", mode: 'copy'

	input:
		tuple val(vamb_key), val(meta), file(fcontigs), path(reads), path(catalogue), path(catalogue_index)

	output:
		tuple val(vamb_key), path(depthout), emit: counttable 
		val(meta), emit: sampleid
		tuple val(meta), file(fcontigs), file(depthout), emit: maps
		tuple val(meta), file(mappingbam), file(mappingbam_index), emit: bam
		path("error.log"),    optional: true, emit: errorlog
		path("versions.yml"), emit: versions

	script:
		sampleID = meta.id

		left_clean = sampleID + "_R1_clean.fastq.gz"
		right_clean = sampleID + "_R2_clean.fastq.gz"
		single_clean = sampleID + "_single_clean.fastq.gz"

		depthout = sampleID + '_depth.txt'
		mappingbam = sampleID + '_mapping_minimap.bam'
		mappingbam_index = sampleID + '_mapping_minimap.bam.bai'

		sample_total_reads = sampleID + '_totalreads.txt'
		if (!meta.single_end) {  
			"""
			#minimap2 -I100G -d catalogue.mmi $catalogue; # make index
			minimap2 -t ${task.cpus} -N 50 -ax sr  $catalogue_index $left_clean $right_clean | samtools view -F 3584 -b --threads ${task.cpus} | samtools sort > $mappingbam 2> error.log # -n 
			samtools index $mappingbam
			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			minimap2: \$(minimap2 --version)
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}' )
			END_VERSIONS

			"""
		} else {
			"""	
			#minimap2 -d catalogue.mmi $catalogue; # make index
			minimap2 -t ${task.cpus} -N 50 -ax sr $catalogue_index $single_clean | samtools view -F 3584 -b --threads ${task.cpus}| samtools sort  > $mappingbam 2> error.log #-n
			samtools index $mappingbam
			jgi_summarize_bam_contig_depths $mappingbam --outputDepth $depthout

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			minimap2: \$(minimap2 --version)
			samtools: \$(samtools --version | head -1 | sed -e "s/samtools //g")
			jgi_summarize_bam_contig_depths: \$(jgi_summarize_bam_contig_depths 2>&1 | head -1 | awk '{print \$2}')
			END_VERSIONS

			"""		
		}
}

process VAMB_COLLECT_DEPTHS {
	cache 'lenient'
	label 'default'
    scratch params.scratch
	tag "$vamb_key"
	//publishDir "${params.outdir}/${sampleID}/vamb", mode: 'copy'

	input:
		tuple val(vamb_key), path(depthout)

	output:
		tuple val(vamb_key), path(alldepths), emit: alldepths
		path("versions.yml"), emit: versions
	script:

		alldepths = vamb_key + '_all_depths.tsv'

		"""
		Rscript ${baseDir}/bin/collectmapping.R $alldepths
		sed -i "s/[.]var/-var/g" $alldepths

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		R: \$(Rscript --version 2>&1 | awk '{print \$4}')
		END_VERSIONS

		"""
}

process VAMB {
	cache 'lenient'
	label 'vamb'
	scratch params.scratch
	tag "$vamb_key"

	input:
        tuple val(vamb_key), path(catalogue), path(catalogue_index), path(alldepths)


	output:
		tuple val(vamb_key), path(cluster_table), emit: all_samples_clustertable
		path("versions.yml"), emit: versions

	script:
		cluster_table = 'all_vamb_contigs_to_bin.tsv'

		"""
		vamb --outdir bin --fasta $catalogue --jgi $alldepths -o _${params.contig_sep}_ 
		mv bin/clusters.tsv $cluster_table

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
		Vamb: 3.0.2
		END_VERSIONS

		"""
}

process VAMB_CONTIGS_SELECTION{
	
	label 'default'
	scratch params.scratch
	tag "$sampleID"
	publishDir "${params.outdir}/vamb/${sampleID}", mode: 'copy', enabled: params.publish_rawbins
	
	input:
		tuple val(vamb_key),val(meta), path(all_cluster_table)
		
	output:
		tuple val(meta), file(persample_clustertable), emit: persample_clustertable

	script:
		sampleID = meta.id

		persample_clustertable = sampleID + '_vamb_contigs_to_bin.tsv'

		"""
		grep ${meta.coassemblygroup}_ $all_cluster_table > $persample_clustertable
		"""
}

process group_vamb {
	
	label 'default'
	scratch params.scratch

    input:
		path(reads_table)
    output:
		path("meta_contigkey.csv"), emit: sample_vambkey
		path("contigs_perkey.csv"), emit: contigs_perkey
		path("temp2_csv.csv"), emit: overview_csv
    script:
    """
    awk '{print int((NR-1)/${params.vamb_groupsize}) "," \$0}' ${reads_table} | sed 's/\\]//' | sed 's/\\[//' > temp2_csv.csv

    #meta and contig-key:
    awk -F, '{print \$2","\$3","\$4","\$1}' temp2_csv.csv > meta_contigkey.csv

    awk -F, '{OFS=","; a[\$1]=a[\$1]" "\$5} END {for (i in a) print i a[i]}' temp2_csv.csv | sed 's/ /,/' > contigs_perkey.csv
    
	"""
}
