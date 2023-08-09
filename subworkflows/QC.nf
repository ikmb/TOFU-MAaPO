/* 
 * Import QC modules
 * MULTIQC1 and MULTIQC2 differ only in having 1 or 2 input-channels, remove this when optional inputs are possible
 */

include {  FILTERREADS } from '../modules/QC/bowtie2'
include {  FASTQC_raw; FASTQC_clean } from '../modules/QC/fastqc'
include {  TRIMREADS; CLEANREADS } from '../modules/QC/bbduk'
include {  COLLECTOR } from '../modules/QC/collect'
/* 
 * QC pipeline logic
 * QC_rmHost removes all Host contaminations, QC_noHost does not filter for Host contaminations
 */
workflow QC_rmHost{
    take: 
        ch_raw_reads
    main:
        if (!params.genome) {
            println "No Host genome was specified!"
        } else if (!params.genomes) {
            exit 1, "Specified a genome name for host mapping, but no genomes are configured for your profile. Aborting!"
        } else if (!params.genomes.containsKey(params.genome)) {
            exit 1, "Specified unknown name for the host genome! Valid presets are: ${params.genomes.keySet()}"
        } else {
            log.info "Using ${params.genome} as host species."

            bowtie_base = params.genomes[params.genome].bowtie_index
        } 
        
        ch_versions = Channel.empty()

        FASTQC_raw(ch_raw_reads)
        TRIMREADS(ch_raw_reads)
        CLEANREADS(TRIMREADS.out.filterReads)
        FILTERREADS(
                CLEANREADS.out.cleanfastq,
                Channel.fromPath("${bowtie_base}*").collect(),         
                Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )  

        ch_cleaned_reads = FILTERREADS.out.cleaned_reads

        FASTQC_clean(ch_cleaned_reads)  

        ch_fastqc_clean_out =  FASTQC_clean.out.fastqc

        ch_versions = ch_versions.mix(TRIMREADS.out.version.first())
        ch_versions = ch_versions.mix(CLEANREADS.out.version.first())
        ch_versions = ch_versions.mix(FILTERREADS.out.version.first())
        ch_versions = ch_versions.mix(FASTQC_raw.out.version.first())
        ch_versions = ch_versions.mix(FASTQC_clean.out.version.first())

    emit:
        qcedreads = ch_cleaned_reads
        fastqcoutput = FASTQC_raw.out.fastqc.collect()
        fastqcoutputclean = ch_fastqc_clean_out.collect()
        versions = ch_versions
}

workflow QC_noHost{
    take: 
        ch_raw_reads
    main:
        
        println "No Host genome was specified, so NO decontamination will be performed!"
        ch_versions = Channel.empty()

        FASTQC_raw(ch_raw_reads)
        TRIMREADS(ch_raw_reads)
        CLEANREADS(TRIMREADS.out.filterReads)
        COLLECTOR(CLEANREADS.out.cleanfastq)  

        ch_cleaned_reads = COLLECTOR.out.cleaned_reads

        FASTQC_clean(ch_cleaned_reads)  

        ch_fastqc_clean_out =  FASTQC_clean.out.fastqc

        ch_versions = ch_versions.mix(TRIMREADS.out.version.first())
        ch_versions = ch_versions.mix(CLEANREADS.out.version.first())
        ch_versions = ch_versions.mix(FASTQC_raw.out.version.first())
        ch_versions = ch_versions.mix(FASTQC_clean.out.version.first())      

    emit:
        qcedreads = ch_cleaned_reads
        fastqcoutput = FASTQC_raw.out.fastqc.collect()
        fastqcoutputclean = ch_fastqc_clean_out.collect()
        versions = ch_versions
}