include { input_check; input_check_qced } from '../subworkflows/input_check'
include { QC_rmHost; QC_noHost } from '../subworkflows/QC'
include { MULTIQC1; MULTIQC2 } from '../subworkflows/QC'
include { metaphlan } from '../subworkflows/metaphlan'
include { kraken } from '../subworkflows/kraken'
include { humann } from '../subworkflows/humann'
include { assembly } from '../subworkflows/assembly'


/* 
 * Main pipeline logic
 */
workflow MW {
    main:

    //QC, either with host contaminations or without:
        if(!params.no_qc){

            input_check()
            ch_raw_reads = input_check.out.reads

            if(params.genome){
                QC_rmHost(
                    ch_raw_reads
                )
                QCout = QC_rmHost.out.qcedreads
                Fastqcoutput = QC_rmHost.out.fastqcoutput.collect()
                FASTqccleanout = QC_rmHost.out.fastqcoutputclean.collect()
            }else{
                QC_noHost(
                    ch_raw_reads
                )
                QCout = QC_noHost.out.qcedreads
                Fastqcoutput = QC_noHost.out.fastqcoutput.collect()
                FASTqccleanout = QC_noHost.out.fastqcoutputclean.collect()
                }
        }else{
            input_check_qced()
            QCout = input_check_qced.out.reads
        }
        
    //kraken:
        if(params.virus || params.kraken || params.bracken){
            kraken(QCout)
        }
    //metaphlan:
        if(params.metaphlan){
            metaphlan(QCout)
        }
    //humann:
        if(params.humann){
            humann(QCout)
        }
    //genome assembly:
        if( params.assembly || params.magscot ){
            assembly(QCout)
        }

    //multiqc, collecting all fastqc- and kraken-files; change this when optional inputs are doable:
        if(!params.no_qc){
            if(params.virus || params.kraken || params.bracken){
                MULTIQC2(
                    Fastqcoutput.collect(),
                    FASTqccleanout.collect(),
                    kraken.out.collect()
                )    
            }else{
                MULTIQC1(
                    Fastqcoutput.collect(),
                    FASTqccleanout.collect()
                ) 
            }
        }
}