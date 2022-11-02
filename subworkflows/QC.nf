/* 
 * Import QC modules
 * MULTIQC1 and MULTIQC2 differ only in having 1 or 2 input-channels, remove this when optional inputs are possible
 */

include {  FILTERREADS_SE;  FILTERREADS_PE } from '../modules/QC/bowtie2'
include {  MULTIQC1;  MULTIQC2 } from '../modules/QC/multiqc'
include {  FASTQC_raw; FASTQC_clean } from '../modules/QC/fastqc'
include {  TRIMREADS; CLEANREADS_PE; CLEANREADS_SE } from '../modules/QC/bbduk'
include {  COLLECTOR_PE; COLLECTOR_SE } from '../modules/QC/collect'
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

     //reads_ch = Channel.fromFilePairs(params.reads , flat: true ).ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
     FASTQC_raw(ch_raw_reads)
     if (!params.single_end) {
         TRIMREADS(ch_raw_reads)
         CLEANREADS_PE(TRIMREADS.out.filterReads)
         CLEANREADS_SE(TRIMREADS.out.filterReads)
         FILTERREADS_PE(
             CLEANREADS_PE.out.join(CLEANREADS_SE.out),
             Channel.fromPath("${bowtie_base}*").collect(),         
             Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )  
         
         ch_cleaned_reads = FILTERREADS_PE.out.cleaned_reads

         FASTQC_clean(ch_cleaned_reads)  

         ch_fastqc_clean_out =  FASTQC_clean.out.fastqc   

     } else {
         TRIMREADS(ch_raw_reads)
         CLEANREADS_SE(TRIMREADS.out.filterReads)
         FILTERREADS_SE(
             CLEANREADS_SE.out,
             Channel.fromPath("${bowtie_base}*").collect(),         
             Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )  
         
         ch_cleaned_reads = FILTERREADS_SE.out.cleaned_reads

         FASTQC_clean(ch_cleaned_reads)  

         ch_fastqc_clean_out =  FASTQC_clean.out.fastqc    
         
     }    
    emit:
     qcedreads = ch_cleaned_reads
     
     fastqcoutput = FASTQC_raw.out.fastqc.collect()
     fastqcoutputclean = ch_fastqc_clean_out.collect()
     
}

workflow QC_noHost{
    take: 
        ch_raw_reads
    main:
        
    println "No Host genome was specified, so NO decontamination will be performed!"

     FASTQC_raw(ch_raw_reads)

     if (!params.single_end) {
         TRIMREADS(ch_raw_reads)
         CLEANREADS_PE(TRIMREADS.out.filterReads)
         CLEANREADS_SE(TRIMREADS.out.filterReads)
         COLLECTOR_PE(CLEANREADS_PE.out.join( CLEANREADS_SE.out) )  
         
         ch_cleaned_reads = COLLECTOR_PE.out.cleaned_reads

         FASTQC_clean(ch_cleaned_reads)  

         ch_fastqc_clean_out =  FASTQC_clean.out.fastqc   

     } else {
         TRIMREADS(ch_raw_reads)
         CLEANREADS_SE(TRIMREADS.out.filterReads)
         COLLECTOR_SE( CLEANREADS_SE.out )  
         
         ch_cleaned_reads = COLLECTOR_SE.out.cleaned_reads

         FASTQC_clean(ch_cleaned_reads)  

         ch_fastqc_clean_out =  FASTQC_clean.out.fastqc    
         
     }         

    emit:
    
     qcedreads = ch_cleaned_reads

     fastqcoutput = FASTQC_raw.out.fastqc.collect()
     fastqcoutputclean = ch_fastqc_clean_out.collect()
}