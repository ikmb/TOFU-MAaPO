process CLEANREADS_PE {
	tag "$sampleID"
	label 'bbmap'

	scratch params.scratch

	input:
		tuple val(meta), path(reads)

	output:
		tuple val(meta), path('*_cleanwithhost.fastq.gz'), emit: cleanfastq
		path('versions.yml'), emit: version
	script:
		sampleID = meta.id
		leftnewname = sampleID + "_1_raw.fastq.gz"
    	rightnewname = sampleID + "_2_raw.fastq.gz"

		left_clean = sampleID + "_R1_cleanwithhost.fastq.gz"
		right_clean = sampleID + "_R2_cleanwithhost.fastq.gz"
		artifact_stats = sampleID + ".bbduk.artifacts.stats"
		
		left_trimmed = sampleID + "_1_trimmed.fastq.gz"
		right_trimmed = sampleID + "_2_trimmed.fastq.gz"
		unpaired = sampleID + "_unpaired_trimmed.fastq.gz"


		"""
		bbduk.sh stats=$artifact_stats threads=${task.cpus} in=${left_trimmed} in2=${right_trimmed} k=31 ref=artifacts,phix ordered cardinality out1=${left_clean} out2=${right_clean} minlength=${params.min_read_length}
		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      		BBMap: \$(bbduk.sh --version 2>&1 | awk 'FNR==2{print \$0}' | sed -e "s/BBMap //g" | sed -e "s/version //g" )
    	END_VERSIONS
		"""
}

process CLEANREADS_SE {
	tag "$sampleID"
	label 'bbmap'

	scratch params.scratch

	input:
		tuple val(meta), path(reads)

	output:
		tuple val(meta),path(unpaired_clean), emit: cleanfastq
		path('versions.yml'), emit: version

	script:
		sampleID = meta.id
		
		left_trimmed = sampleID + "_1_trimmed.fastq.gz"
		right_trimmed = sampleID + "_2_trimmed.fastq.gz"
		unpaired = sampleID + "_unpaired_trimmed.fastq.gz"

		unpaired_clean = sampleID + "_single_cleanwithhost.fastq.gz"

		"""
		bbduk.sh threads=${task.cpus} in=$unpaired  k=31 ref=artifacts,phix ordered cardinality out1=${unpaired_clean} minlength=${params.min_read_length}
		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	  BBMap: \$(bbduk.sh --version 2>&1 | awk 'FNR==2{print \$0}' | sed -e "s/BBMap //g" | sed -e "s/version //g" )
    	END_VERSIONS
		"""
}

process TRIMREADS {
	tag "$sampleID"
	label 'bbmap'
	//errorStrategy 'ignore'
	scratch params.scratch

	input:
		tuple val(meta), path(reads)
		
	output:
		
		tuple val(meta), path('*_trimmed.fastq.gz'), emit: filterReads
		path('versions.yml'), emit: version
	script:
		sampleID = meta.id

		bbduk_adapter_stats = sampleID + ".bbduk.adapter.stats"

		leftnewname = sampleID + "_1_raw.fastq.gz"
    	rightnewname = sampleID + "_2_raw.fastq.gz"

		left_trimmed = sampleID + "_1_trimmed.fastq.gz"
		right_trimmed = sampleID + "_2_trimmed.fastq.gz"
		unpaired = sampleID + "_unpaired_trimmed.fastq.gz"

		if (meta.single_end) {
		"""
		[ ! -f  $leftnewname ] && ln -s ${reads} $leftnewname

		bbduk.sh stats=$bbduk_adapter_stats \
				threads=${task.cpus} \
				in=${leftnewname} \
				out=${left_trimmed} \
				outs=$unpaired \
				ref=${params.adapters} \
				ktrim=r \
				k=23 \
				mink=11 \
				hdist=1 \
				minlength=${params.min_read_length} \
				tpe \
				tbo
		rm ${left_trimmed}

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	  BBMap: \$(bbduk.sh --version 2>&1 | awk 'FNR==2{print \$0}' | sed -e "s/BBMap //g" | sed -e "s/version //g" )
    	END_VERSIONS
		"""
		} else {
		"""
		[ ! -f  $leftnewname ] && ln -s ${reads[0]} $leftnewname
    	[ ! -f  $rightnewname ] && ln -s ${reads[1]} $rightnewname

		bbduk.sh stats=$bbduk_adapter_stats \
				threads=${task.cpus} \
				in=${leftnewname} \
				in2=${rightnewname} \
				out1=${left_trimmed} \
				out2=${right_trimmed} \
				outs=$unpaired \
				ref=${params.adapters} \
				ktrim=r \
				k=23 \
				mink=11 \
				hdist=1 \
				minlength=${params.min_read_length} \
				tpe \
				tbo

		cat <<-END_VERSIONS > versions.yml
    	"${task.process}":
      	  BBMap: \$(bbduk.sh --version 2>&1 | awk 'FNR==2{print \$0}' | sed -e "s/BBMap //g" | sed -e "s/version //g" )
    	END_VERSIONS
		"""
		}
}