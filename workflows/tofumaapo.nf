include { input_check; input_check_qced; input_sra } from '../subworkflows/input_check'
include { QC } from '../subworkflows/QC'
include { MULTIQC } from '../modules/QC/multiqc'
include { metaphlan } from '../subworkflows/metaphlan'
include { kraken } from '../subworkflows/kraken'
include { humann } from '../subworkflows/humann'
include { salmon } from '../subworkflows/salmon'
include { assembly } from '../subworkflows/assembly'
include { SOFTWARE_VERSIONS } from '../modules/software_versions'


/* 
 * Main pipeline logic
 */

workflow tofumaapo {
	main:

		ch_versions = Channel.from([])

	// inputs:
		if(!params.no_qc){
			if(params.reads && params.sra){
				exit 1, "Please only declare either --sra or --read. Not both!"
			}

			if(params.reads){
			input_check()
			ch_raw_reads = input_check.out.reads
			} else {
				if(params.sra){
					input_sra()
					ch_raw_reads = input_sra.out.reads
				}else{
					exit 1, "No input was declared! Please declare input with --reads or --sra !"
				}
			}
	//QC:
			QC(ch_raw_reads)
			QCout = QC.out.qcedreads
			Fastqcoutput = QC.out.qcreports
			ch_versions = ch_versions.mix( QC.out.versions )

		}else{
			input_check_qced()
			QCout = input_check_qced.out.reads
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

	//metaphlan:
		if(params.metaphlan || params.updatemetaphlan){
			metaphlan(QCout)

			ch_versions = ch_versions.mix( metaphlan.out.versions )
		}

	//humann:
		if(params.humann || params.updatehumann){
			humann(QCout)

			ch_versions = ch_versions.mix( humann.out.versions )
		}
		
	//genome assembly:
		if( params.assembly || params.coassembly ){
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