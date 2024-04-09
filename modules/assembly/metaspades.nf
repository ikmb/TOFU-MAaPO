process METASPADES {

	publishDir {"${params.outdir}/${meta.id}/metaspades"}, mode: 'copy', pattern: '**/*final.contigs.fa', enabled: params.publish_contigs
	scratch params.scratch
	label 'spades'
	tag "$meta.id"

	input:
		tuple val(meta), path(reads)

	output:
		path("**/*"), emit: outputfolder
		tuple val(coassemblygroup), file(output_final_contigs), emit: contigs
		//tuple val(meta), file(output_final_contigs), path('*_clean.fastq.gz', includeInputs: true), emit: contigs
		path("versions.yml"),          optional: true, emit: versions
	script:
		def memory = task.memory.toGiga()
		def input = meta.single_end ? "--12  ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]}"
		coassemblygroup = meta.coassemblygroup
		replacer = ">" + coassemblygroup + "_k"
		output_final_contigs = coassemblygroup + "_spades_final.contigs.fa"
		contig_sep = 'metaspadescontig'
		"""
		spades.py --meta $input -k 21,33,55 -o spades_out -t ${task.cpus} -m ${memory}
		#convert the contigs.fasta file into the layout of megahit outputs. Give new names for contigs, change order of contig data, so that length is last and remove newline characters within the contigs sequences
		awk '{gsub (/^>NODE_/, "$replacer");print}' spades_out/contigs.fasta | awk '{OFS=""} {gsub(/_k/, "_${contig_sep}_");gsub(/_cov_/, " cov=");gsub(/_length_/, " len=");print}' | awk '{FS=" "} (\$1 ~ />/) {print \$1" "\$3" placeholder=1 "\$2}  (\$1 !~ />/){print \$0}' | awk '!/^>/ {printf "%s", \$0; next} {if (!first) {first=1} else {print ""} print}' > $output_final_contigs


		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		Python: \$(python3 --version | sed -e "s/Python //g" )
		END_VERSIONS

		"""
	}

/*
zcat $unpaired > unpaired.fq
zcat $left > left.fq
zcat $right > right.fq
*/