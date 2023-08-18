process MEGAHIT {
	scratch params.scratch
	label 'megahit'
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

process MEGAHIT_coassembly {
	//scratch params.scratch
	scratch false
	label 'megahit'
	//tag "$sampleID"
	//publishDir {"${params.outdir}/Megahit"}, mode: 'copy', pattern: '*_final.contigs.fa', enabled: params.publish_megahit

	input:
		tuple val(meta), path(reads)

	output:
		path('output/*'), emit: outputfolder
		//tuple val(meta), file(output_final_contigs), path('*_clean.fastq.gz', includeInputs: true), emit: contigs
		tuple val(meta), file(output_final_contigs), emit: contigs
        path("versions.yml"),          optional: true, emit: versions

	shell:
		//sampleID = meta.id

		replacer = ">" + meta + "_k"
		output_final_contigs = meta + "_final.contigs.fa"
        //def input = "-r \"" + single.join(",") + "\"" + "-1 \"" + links.join(",") + "\" -2 \"" + rechts.join(",") + "\""
		
		"""
        #!/bin/bash

        file_paths="${reads.flatten().join(",")}"

        r1_files=""
        r2_files=""
        unpaired_files=""

        IFS=',' read -ra files <<< "\$file_paths"

        for file in "\${files[@]}"; do
            if [[ \$file == *_R1* ]]; then
                r1_files+="\$file,"
            elif [[ \$file == *_R2* ]]; then
                r2_files+="\$file,"
            else
                unpaired_files+="\$file,"
            fi
        done

        r1_files=\${r1_files%,}  # Remove trailing comma
        r2_files=\${r2_files%,}  # Remove trailing comma
        unpaired_files=\${unpaired_files%,}  # Remove trailing comma

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

		awk '{gsub (/^>k/, "$replacer");print}' output/final.contigs.fa | awk '{OFS=""} {gsub(/_k[0-9]+_/, "_${params.contig_sep}_");print}' > $output_final_contigs

		cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        MEGAHIT: \$(megahit --version 2>&1 | sed -e "s/MEGAHIT v//g")
        END_VERSIONS

        """	
}
//-m 0.95 \