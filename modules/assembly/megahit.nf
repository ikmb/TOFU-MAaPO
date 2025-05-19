//Deprecated, only for single sample assembly:
process MEGAHIT {
	scratch params.scratch
	label 'megahit'
	label 'very_long_run'
	tag "$sampleID"
	publishDir {"${params.outdir}/Megahit"}, mode: 'copy', pattern: '*_final.contigs.fa', enabled: params.publish_megahit

	input:
		tuple val(meta), path(reads)

	output:
		path('output/*'), emit: outputfolder
		tuple val(meta), file(output_final_contigs), path('*_clean.fastq.gz', includeInputs: true), emit: contigs
		path("versions.yml"),          optional: true, emit: versions

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

		if (!meta.single_end) {  
			"""
			zcat $unpaired_clean > $fq_single
			zcat $left_clean > $fq_left
			zcat $right_clean > $fq_right
			
			megahit -1 $fq_left \
				-2 $fq_right \
				-r $fq_single \
				-o output \
				-m ${task.memory.toBytes()} \
				-t ${task.cpus}	

			rm $fq_single
			rm $fq_left
			rm $fq_right

			awk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa | awk '{OFS=""} {gsub(/_k[0-9]+_/, "_${params.contig_sep}_");print}' > $output_final_contigs

			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			MEGAHIT: \$(megahit --version 2>&1 | sed -e "s/MEGAHIT v//g")
			END_VERSIONS
		
			"""
		} else {
			"""	
			zcat $unpaired_clean > $fq_single

			megahit	-r $fq_single \
				-o output \
				-m ${task.memory.toBytes()} \
				-t ${task.cpus}	

			rm $fq_single

			#awk '{gsub (/^>k/, "/$replacer");print}' output/final.contigs.fa > $output_final_contigs
			awk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa | awk '{OFS=""} {gsub(/_k[0-9]+_/, "_${params.contig_sep}_");print}' > $output_final_contigs
			
			cat <<-END_VERSIONS > versions.yml
			"${task.process}":
			MEGAHIT: \$(megahit --version 2>&1 | sed -e "s/MEGAHIT v//g")
			END_VERSIONS
			
			"""			
		}
}
//-m 0.95 \

//Robust version that sorts reads by R1, R2 or unpaired. Used for coassembly and single sample assembly
process MEGAHIT_assembly {
	scratch params.scratch
	label 'megahit'
	label 'very_long_run'
	tag "$coassemblygroup"
	publishDir {"${params.outdir}/Megahit"}, mode: 'copy', pattern: '*_final.contigs.fa', enabled: params.publish_megahit

	input:
		tuple val(coassemblygroup), path(reads)

	output:
		path('output/*'), emit: outputfolder
		tuple val(coassemblygroup), file(output_final_contigs), emit: contigs
		path("versions.yml"),          optional: true, emit: versions

	shell:
		coassemblygroup
		
		replacer = ">" + coassemblygroup + "_k"
		output_final_contigs = coassemblygroup + "_final.contigs.fa"
		//def input = "-r \"" + single.join(",") + "\"" + "-1 \"" + links.join(",") + "\" -2 \"" + rechts.join(",") + "\""
		
		"""
		#!/bin/bash

		file_paths="${reads.flatten().join(",")}"

		#sort read files based on their _R[1,2]_ flag
		r1_files=""
		r2_files=""
		unpaired_files=""

		IFS=',' read -ra files <<< "\$file_paths"

		for file in "\${files[@]}"; do
			if [[ \$file == *_R1_* ]]; then
				r1_files+="\$file,"
			elif [[ \$file == *_R2_* ]]; then
				r2_files+="\$file,"
			else
				unpaired_files+="\$file,"
			fi
		done

		# Remove trailing comma
		r1_files=\${r1_files%,}  
		r2_files=\${r2_files%,}
		unpaired_files=\${unpaired_files%,}

		echo "R1 Files: \$r1_files"
		echo "R2 Files: \$r2_files"

		if [[ -n "\$unpaired_files" ]]; then
			echo "Unpaired Files: \$unpaired_files"

			megahit -1 \$r1_files \
				-2 \$r2_files \
				-r \$unpaired_files \
				-o output \
				-m ${task.memory.toBytes()} \
				-t ${task.cpus}	

		else
			echo "No other files."
			megahit	-1 \$r1_files \
				-2 \$r2_files \
				-o output \
				-m ${task.memory.toBytes()} \
				-t ${task.cpus}	
		fi

		[ -s output/final.contigs.fa ] || { echo "ERROR - final.contigs.fa is empty. Exiting."; exit 1; }

		awk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa | awk '{OFS=""} {gsub(/_k[0-9]+_/, "_${params.contig_sep}_");print}' > $output_final_contigs

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
		MEGAHIT: \$(megahit --version 2>&1 | sed -e "s/MEGAHIT v//g")
		END_VERSIONS

		"""	
}
