process HUMANN {
    label 'humann'
    scratch true
    //publishDir "${params.outdir}/${sampleID}/humann", mode: 'copy'

     input:
     tuple val(sampleID),path(left_clean),path(right_clean),path(unpaired_clean)
     output:
	 path(genefamilies), emit: genefamilies
     path(pathabundance), emit: pathabundance
     path(pathcoverage), emit: pathcoverage
     script:
     genefamilies = sampleID + '_genefamilies.tsv'
     pathabundance = sampleID + '_pathabundance.tsv'
     pathcoverage = sampleID + '_pathcoverage.tsv'
	 merged = sampleID + '.fq'
     """
		zcat $left_clean > left.fq
    	zcat $right_clean > right.fq
    	cat left.fq right.fq > $merged
    	humann --input $merged --output . --remove-temp-output --threads ${task.cpus} --nucleotide-database ${params.humann_db}/chocophlan --protein-database ${params.humann_db}/uniref --metaphlan-options "--bowtie2db ${params.metaphlan_db} -x mpa_v30_CHOCOPhlAn_201901 --nproc ${task.cpus}"
        rm *.fq 
     """
	}
    //tuple path(genefamilies),path(pathabundance),path(pathcoverage)
    //tuple val(sampleID),file("${sampleID}_genefamilies.tsv"),file("${sampleID}_pathabundance.tsv"),file("${sampleID}_pathcoverage.tsv")
	process JOINgenefamilies {
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    label 'humann'

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
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    label 'humann'

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
    publishDir "${params.outdir}/humann", mode: 'copy'
    scratch true
    label 'humann'

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