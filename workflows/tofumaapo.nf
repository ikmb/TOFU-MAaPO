include { input_check; input_check_qced; input_sra } from '../subworkflows/input_check'
include { QC_rmHost; QC_noHost } from '../subworkflows/QC'
include { MULTIQC1; MULTIQC2 } from '../subworkflows/QC'
include { metaphlan } from '../subworkflows/metaphlan'
include { kraken } from '../subworkflows/kraken'
include { humann } from '../subworkflows/humann'
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
                    exit 1, "No input in --reads or --sra was declared!"
                    }
            }
    //QC, either with host contaminations or without:
            if(params.genome){
                QC_rmHost(
                    ch_raw_reads
                )
                QCout = QC_rmHost.out.qcedreads
                Fastqcoutput = QC_rmHost.out.fastqcoutput.collect()
                FASTqccleanout = QC_rmHost.out.fastqcoutputclean.collect()

                ch_versions = ch_versions.mix( QC_rmHost.out.versions )
            }else{
                QC_noHost(
                    ch_raw_reads
                )
                QCout = QC_noHost.out.qcedreads
                Fastqcoutput = QC_noHost.out.fastqcoutput.collect()
                FASTqccleanout = QC_noHost.out.fastqcoutputclean.collect()

                ch_versions = ch_versions.mix( QC_noHost.out.versions )
                }
        }else{
            input_check_qced()
            QCout = input_check_qced.out.reads
        }
        
    //kraken:
        if(params.virus || params.kraken || params.bracken){
            kraken(QCout)

            ch_versions = ch_versions.mix( kraken.out.versions )
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
        if( params.assembly || params.magscot ){
            assembly(QCout)

            ch_versions = ch_versions.mix( assembly.out.versions )
        }

    //multiqc, collecting all fastqc- and kraken-files; change this when optional inputs are doable:
        if(!params.no_qc){
            if(params.virus || params.kraken || params.bracken){
                MULTIQC2(
                    Fastqcoutput.collect(),
                    FASTqccleanout.collect(),
                    kraken.out.kraken_data.collect()
                )    
            }else{
                MULTIQC1(
                    Fastqcoutput.collect(),
                    FASTqccleanout.collect()
                ) 
            }
        }

    SOFTWARE_VERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}