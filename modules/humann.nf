
    process HUMANN {

    label 'humann'
    tag "$sampleID"
    scratch params.scratch
    //scratch true
    //publishDir "${params.outdir}/${sampleID}/humann", mode: 'copy'

     input:
     tuple val(meta),path(reads)

     output:
	 path('*_genefamilies.tsv'), emit: genefamilies
     path('*_pathabundance.tsv'), emit: pathabundance
     path('*_pathcoverage.tsv'), emit: pathcoverage
     path('*'), emit: humannouts

     script:
     sampleID = meta.id

	left_clean = sampleID + "_R1_clean.fastq.gz"
	right_clean = sampleID + "_R2_clean.fastq.gz"
	unpaired_clean = sampleID + "_single_clean.fastq.gz"

    phlan_left = sampleID + "_1.fq"
    phlan_right = sampleID + "_2.fq"
    phlan_single = sampleID + "_single.fq"

    genefamilies = sampleID + '_genefamilies.tsv'
    pathabundance = sampleID + '_pathabundance.tsv'
    pathcoverage = sampleID + '_pathcoverage.tsv'
	merged = sampleID + '.fq'

    if (!params.single_end) {
        """
		zcat ${left_clean} > $phlan_left
    	zcat ${right_clean} > $phlan_right
        
        #check if unpaired/single reads are present
        if [ -s ${unpaired_clean} ];then

            zcat ${unpaired_clean} > $phlan_single
            cat $phlan_left $phlan_right $phlan_single > $merged
        else
            cat $phlan_left $phlan_right > $merged
        fi

    	humann --input $merged \
            --output . \
            --remove-temp-output \
            --threads ${task.cpus} \
            --nucleotide-database ${params.humann_db}/chocophlan \
            --protein-database ${params.humann_db}/uniref \
            --metaphlan-options "--bowtie2db ${params.metaphlan_db} \
            -x mpa_v30_CHOCOPhlAn_201901 \
            --stat_q 0.2 \
            --force \
            -t rel_ab_w_read_stats \
            --nproc ${task.cpus}" \
            --output-basename $sampleID
        rm *.fq 
        """
	} else {
        """
        zcat ${unpaired_clean} > $phlan_single

    	humann --input $phlan_single \
            --output . \
            --output-basename $sampleID \
            --remove-temp-output \
            --threads ${task.cpus} \
            --nucleotide-database ${params.humann_db}/chocophlan \
            --protein-database ${params.humann_db}/uniref \
            --metaphlan-options "--bowtie2db ${params.metaphlan_db} \
            -x mpa_v30_CHOCOPhlAn_201901 \
            --stat_q 0.2 \
            --force \
            -t rel_ab_w_read_stats \
            --nproc ${task.cpus}"
        rm *.fq 
        """
    }
}

    //tuple path(genefamilies),path(pathabundance),path(pathcoverage)
    //tuple val(sampleID),file("${sampleID}_genefamilies.tsv"),file("${sampleID}_pathabundance.tsv"),file("${sampleID}_pathcoverage.tsv")
	process JOINgenefamilies {
    
    label 'humann'
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    
    input:
	path('*')

    output:
        file(mergedtable)

    script:
        mergedtable = "humann_merged_genefamilies.tsv"

        """
            humann_join_tables --input . --output $mergedtable
        """
  	}

	process JOINpathabundance {

    label 'humann'
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    

    input:
	path('*')

    output:
        file(mergedtable)

    script:
        mergedtable = "humann_merged_pathabundance.tsv"

        """
            humann_join_tables --input . --output $mergedtable
        """
  	}

	process JOINpathcoverage {

    label 'humann'
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    

    input:
	path('*')

    output:
        file(mergedtable)

    script:
        mergedtable = "humann_merged_pathcoverage.tsv"

        """
            humann_join_tables --input . --output $mergedtable
        """
  	}       