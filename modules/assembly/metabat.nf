process METABAT {

    label 'metabat2'
    scratch params.scratch
    tag "${sampleID}_${meta.assembler}"
    publishDir "${params.outdir}/Metabat2/${sampleID}", mode: 'copy', enabled: params.publish_rawbins

    input: 
        tuple val(meta), file(fcontigs), file(depthout)

    output:
		tuple val(meta), file("${sampleID}_bin.*.fa"), optional: true, emit:metabatout
		path("versions.yml"),          optional: true, emit: versions

	script:
		sampleID = meta.id
		"""
		metabat2 -i $fcontigs -a ${depthout} -o ${sampleID}_bin -t ${task.cpus} -m ${params.contigsminlength}

		touch ${sampleID}_bin.emptydummy.fa

		cat <<-END_VERSIONS> versions.yml
		"${task.process}":
		METABAT: \$(metabat2 --help 2>&1 | awk 'NR==2{print}' | sed -n -e 's/^.*version //p' | awk '{print \$1}')
		END_VERSIONS
		
		"""
	stub:
		sampleID = meta.id
		"""
		touch ${sampleID}_bin.stub.fa

		echo "METABAT_stub" > versions.yml
		"""
}

process contigs_to_bins {

	label 'default'
	scratch params.scratch
	tag "${sampleID}_${meta.assembler}"
	publishDir {"${params.outdir}/Metabat2/${sampleID}"}, mode: 'copy', enabled: params.publish_rawbins as boolean
	
	input: 
		tuple val(meta), file(fafile)

	output:
		tuple val(meta), file("${sampleID}_metabat2_magscot_contigs_to_bin.tsv"), emit: metabat2_contigs_to_bins

	script:
		sampleID = meta.id
		
		if (fafile instanceof List && fafile.size() > 1) {
			"""
			echo "Variable contains multiple files"
			grep '>' ${fafile} | tr '>' '\t' | tr -d ':' > ${sampleID}_metabat2_contigs_to_bin.tsv
			
			gawk '{print \$1"\t"\$2"\tmetabat2"}'  ${sampleID}_metabat2_contigs_to_bin.tsv > ${sampleID}_metabat2_magscot_contigs_to_bin.tsv
			"""
		} else {
			"""
			echo "Variable contains a single file"
			set +e
			{
				grep '>' ${fafile} | tr '>' '\t' | tr -d ':' | awk -F'\t' '{OFS="\t"; if (\$1=="") \$1="${fafile}"; print \$0 }' > ${sampleID}_metabat2_contigs_to_bin.tsv
			}
			set -e

			gawk '{print \$1"\t"\$2"\tmetabat2"}'  ${sampleID}_metabat2_contigs_to_bin.tsv > ${sampleID}_metabat2_magscot_contigs_to_bin.tsv
			"""
		}
	stub:
		sampleID = meta.id
		"""
		touch ${sampleID}_metabat2_magscot_contigs_to_bin.tsv

		echo "contigs_to_bins_stub" > versions.yml
		"""
}