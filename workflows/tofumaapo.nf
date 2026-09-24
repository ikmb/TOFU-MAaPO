include { input_check; input_sra } from '../subworkflows/input_check'
include { QC } from '../subworkflows/QC'
include { MIN_READS as MIN_RAW_READS } from '../modules/QC/min_reads'
include { MIN_READS as MIN_QCED_READS } from '../modules/QC/min_reads'
include { MULTIQC } from '../modules/QC/multiqc'
include { metaphlan } from '../subworkflows/metaphlan'
include { kraken } from '../subworkflows/kraken'
include { humann } from '../subworkflows/humann'
include { salmon } from '../subworkflows/salmon'
include { sylph } from '../subworkflows/sylph'
include { assembly } from '../subworkflows/assembly'
include { SOFTWARE_VERSIONS } from '../modules/software_versions'


/* 
 * Main pipeline logic
 */

workflow tofumaapo {
	main:

		ch_versions = Channel.from([])
		ch_input_reads = Channel.from([])

	// inputs:
		
		if(!params.reads && !params.sra){
			exit 1, "Please declare an input! You can do so with --sra and --reads."
		}

		if(params.reads){
			input_check()
			ch_input_reads = ch_input_reads.mix(input_check.out.reads)
		}

		if(params.sra){
			input_sra()
			ch_input_reads = ch_input_reads.mix(input_sra.out.reads)
		}

	//QC:
		if(params.min_raw_reads > 0) {
			MIN_RAW_READS(ch_input_reads, params.min_raw_reads, 'raw')
			ch_raw_reads = MIN_RAW_READS.out.fastq
		}else{
			ch_raw_reads = ch_input_reads
		}
		
		if(!params.no_qc){			
			QC(ch_raw_reads)
			Fastqcoutput = QC.out.qcreports
			ch_versions = ch_versions.mix( QC.out.versions )

			if(params.min_qced_reads > 0) {
				MIN_QCED_READS(QC.out.qcedreads, params.min_qced_reads, 'QCed')
				QCout = MIN_QCED_READS.out.fastq
			}else{
				QCout = QC.out.qcedreads
			}
		}else{
			QCout = ch_raw_reads
		}

	//kraken:
		if(params.virus || params.kraken || params.bracken){
			kraken(QCout)

			ch_versions = ch_versions.mix( kraken.out.versions )
		}

	//salmon:
		if(params.salmon){
			salmon(QCout)

			ch_versions = ch_versions.mix( salmon.out.versions)
		}

	//sylph
		if(params.sylph){
			sylph(QCout)

			ch_versions = ch_versions.mix (sylph.out.versions)
		}

	//metaphlan:
		if(params.metaphlan || params.updatemetaphlan){
			metaphlan(QCout)

			ch_versions = ch_versions.mix( metaphlan.out.versions )
		}

	//humann:
		if(params.humann || params.updatehumann){
			if(params.updatemetaphlan){
				mphlan_ready = metaphlan.out.metaphlan_ready
			}else{
				mphlan_ready = Channel.of('true')
			}
			humann(QCout, mphlan_ready)

			ch_versions = ch_versions.mix( humann.out.versions )
		}
		
	//genome assembly:
		if( params.assembly || params.coassembly || params.updategtdbtk ){
			assembly(QCout)

			ch_versions = ch_versions.mix( assembly.out.versions )
		}

	//multiqc, collecting all fastqc- and kraken-files
		if(!params.no_qc){
			if(params.virus || params.kraken || params.bracken){
				ch_multiqcinputs = Fastqcoutput.collect().mix( kraken.out.kraken_data.collect() ).collect()
			}else{
				ch_multiqcinputs = Fastqcoutput.collect()
			}

			MULTIQC( ch_multiqcinputs )
			ch_versions = ch_versions.mix( MULTIQC.out.versions )
		}

	SOFTWARE_VERSIONS (
		ch_versions.unique().collectFile(name: 'collated_versions.yml')
	)
}