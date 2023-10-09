/* 
 * Import QC modules
 */

include {  FILTERREADS } from '../modules/QC/bowtie2'
include {  FASTQC_raw; FASTQC_clean } from '../modules/QC/fastqc'
include {  TRIMREADS; CLEANREADS } from '../modules/QC/bbduk'
include {  COLLECTOR } from '../modules/QC/collect'
include {  FASTP } from '../modules/QC/fastp'

/* 
 * QC pipeline logic
 */

workflow QC{
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

		if(!params.fastp){
			FASTQC_raw(ch_raw_reads)
			TRIMREADS(ch_raw_reads)
			CLEANREADS(TRIMREADS.out.filterReads)

			ch_versions = ch_versions.mix(FASTQC_raw.out.version.first())
			ch_versions = ch_versions.mix(TRIMREADS.out.version.first())
			ch_versions = ch_versions.mix(CLEANREADS.out.version.first())
			if(!params.genome){
				COLLECTOR(CLEANREADS.out.cleanfastq) 

				ch_cleaned_reads = COLLECTOR.out.cleaned_reads
			}else{
				FILTERREADS(
					CLEANREADS.out.cleanfastq,
					Channel.fromPath("${bowtie_base}*").collect(),         
					Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )  

				ch_versions = ch_versions.mix(FILTERREADS.out.version.first())

				ch_cleaned_reads = FILTERREADS.out.cleaned_reads
			}

			FASTQC_clean(ch_cleaned_reads)
			ch_versions = ch_versions.mix(FASTQC_clean.out.version.first())

			ch_fastqc_clean_out =  FASTQC_clean.out.fastqc

			ch_qcreports = FASTQC_raw.out.fastqc.collect().mix( ch_fastqc_clean_out.collect() )
		}else{
			FASTP(ch_raw_reads)
			ch_versions = ch_versions.mix(FASTP.out.version.first())

			ch_qcreports = FASTP.out.html_report.collect()
			
			if(!params.genome){
				COLLECTOR(FASTP.out.filterReads) 

				ch_cleaned_reads = COLLECTOR.out.cleaned_reads
			}else{
				FILTERREADS(
					FASTP.out.filterReads,
					Channel.fromPath("${bowtie_base}*").collect(),
					Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )  

				ch_versions = ch_versions.mix(FILTERREADS.out.version.first())

				ch_cleaned_reads = FILTERREADS.out.cleaned_reads

				FASTQC_clean(ch_cleaned_reads)
				ch_versions = ch_versions.mix(FASTQC_clean.out.version.first())
			}
		}
		
	emit:
		qcedreads = ch_cleaned_reads
		qcreports = ch_qcreports
		versions = ch_versions
}