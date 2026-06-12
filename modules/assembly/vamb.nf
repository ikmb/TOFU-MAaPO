process VAMB_CATALOGUE {
	label 'vamb'
	label 'long_run'
	scratch params.scratch
	tag "$vamb_key"

	input:
		tuple val(vamb_key), path(contigs)

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

process VAMB_strobealign {
	label 'strobealing'
	label 'long_run'
	scratch params.scratch
	tag "${vamb_key}_${meta.id}"
	
	input:
		tuple val(vamb_key), val(meta), path(reads), path(catalogue)

	output:
		tuple val(vamb_key), path(abundance_table), emit: abundance
		path("versions.yml"), emit: versions
	script:
		sampleID = meta.id
		abundance_table = sampleID + "_abundance.tsv"
		//strobealign cannot deal with the unpaired file in a paired read set
		def selected_reads = reads.flatten().collect().size() >= 2 ? reads[0..1] : [reads[0]]
    	def reads_str = selected_reads.collect { read -> read.toString() }.join(' ')
		"""
		strobealign -t ${task.cpus} --aemb $catalogue $reads_str > ${abundance_table}
		
		cat <<-END_VERSIONS> versions.yml
		"${task.process}":
		strobealign: \$(strobealign --version )
		END_VERSIONS
		"""
}

process VAMB_merge_aemb {
	label 'vamb'
	label 'short_run'
	scratch params.scratch
	tag "$vamb_key"

	input:
		tuple val(vamb_key), path(abundance_table)
	output:
		tuple val(vamb_key), path(combined_abundance), emit: abundance
		path("versions.yml"), emit: versions
	script:
		combined_abundance = vamb_key + "_abundance.tsv"

		"""
		mkdir input && mv ${abundance_table} input/

		python3 /workspace/vamb/src/merge_aemb.py input/ $combined_abundance

		cat <<-END_VERSIONS> versions.yml
		"${task.process}":
		Python: \$(python --version | sed -e "s/Python //g" )
		END_VERSIONS

		"""
	stub:
		alldepths = vamb_key + '_all_depths.tsv'
		"""
		touch $alldepths
		echo "VAMB_COLLECT_DEPTHS_stub" > versions.yml
		"""
}

process VAMB {
	cache 'lenient'
	label 'vamb'
	label 'gpu'
	label 'exclusive' // vamb does not control for numpy threads, which takes all threads by default
	scratch params.scratch
	tag "$vamb_key"
	containerOptions { params.gpu ? '--nv' : '' }
	input:
        tuple val(vamb_key), path(catalogue), path(abundance)//path(alldepths)

	output:
		tuple val(vamb_key), path(cluster_table), emit: all_samples_clustertable
		path("versions.yml"), emit: versions

	script:
		cluster_table = 'all_vamb_contigs_to_bin.tsv'
		def gpu_option = params.gpu ? "--cuda" : ""
		"""
		vamb bin default \
			--outdir bin \
			--fasta $catalogue \
			--abundance_tsv $abundance \
			-p ${task.cpus} \
			$gpu_option \
			-o _${params.contig_sep}_
			
		cat bin/vae_clusters_split.tsv > $cluster_table

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
      	Python: \$(python --version | sed -e "s/Python //g" )
		Vamb: 3.0.2
		END_VERSIONS

		"""
	stub:
		cluster_table = 'all_vamb_contigs_to_bin.tsv'
		"""
		touch $cluster_table
		echo "VAMB_stub" > versions.yml
		"""
}

process VAMB_CONTIGS_SELECTION{
	label 'short_run'
	label 'default'
	scratch params.scratch
	tag "$sampleID"
	//publishDir "${params.outdir}/vamb/${sampleID}", mode: 'copy', enabled: params.publish_rawbins
	publishDir "${params.outdir}/vamb", mode: 'copy', enabled: params.publish_rawbins,
        saveAs: { filename -> "${meta.id}/${filename}" }
	
	input:
		tuple val(vamb_key),val(meta), path(all_cluster_table)
		
	output:
		tuple val(meta), file(persample_clustertable), optional: true, emit: persample_clustertable
		tuple val(meta), file(formatted_contigs_to_bin),optional: true, emit: magscot_contigbinlist
	script:
		sampleID = meta.id

		persample_clustertable = sampleID + '_vamb_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_vamb_magscot_contigs_to_bin.tsv'
		"""
		grep ${meta.coassemblygroup}_ $all_cluster_table > $persample_clustertable

		gawk '{print \$1"\t"\$2"\tvamb"}'  $persample_clustertable > $formatted_contigs_to_bin
		"""
	stub:
		sampleID = meta.id
		persample_clustertable = sampleID + '_vamb_contigs_to_bin.tsv'
		formatted_contigs_to_bin = sampleID + '_vamb_magscot_contigs_to_bin.tsv'
		"""
		touch $persample_clustertable
		touch $formatted_contigs_to_bin
		echo "VAMB_CONTIGS_SELECTION_stub" > versions.yml
		"""
}

process group_vamb {
	label 'default'
	label 'short_run'
	scratch params.scratch
	publishDir "${params.outdir}/", mode: 'copy', enabled: params.publish_rawbins
	
	input:
		path(reads_table)
	output:
		path("binninggroup_to_sample.csv"), emit: sample_vambkey
		//path("contigs_perkey.csv"), emit: contigs_perkey
		//path("temp2.csv"), emit: overview_csv
	script:
	"""
		# If all values for co-binning are unique, create new VAMB groups.
		if awk  -F';' 'NF { if (++seen[\$1] > 1) dup=1 } END { exit dup }' "${reads_table}"; then
			echo "unique grouping"
			awk -F';' '{print int((NR-1)/${params.vamb_groupsize}) ";" \$2}' ${reads_table} > binninggroup_to_sample.csv
		else
			# Co-binning is to be performed with the given co-binning grouping.
			echo "using available grouping"
			awk -F';' '{print \$1 ";" \$2}' ${reads_table} > binninggroup_to_sample.csv
		fi
		"""
	}
//#awk -F';' '{a[\$1]=a[\$1]"\t"\$4} END {for (i in a) print i a[i]}' temp2.csv | sed 's/ /,/' > contigs_perkey.csv
	
