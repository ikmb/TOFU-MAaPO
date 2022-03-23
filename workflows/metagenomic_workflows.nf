include { input_check } from '../subworkflows/input_check'
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
        input_check()
        ch_raw_reads = input_check.out.reads
    //QC, either with host contaminations or without:
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
//        Channel.empty().join(QC_rmHost.out.qcedreads).view()
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
        if(params.assembly){
            assembly(QCout)
        }

    //multiqc, collecting all fastqc- and kraken-files; change this when optional inputs are doable:
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