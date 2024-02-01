include { input_check; input_sra; input_sra_list } from '../subworkflows/input_check'
include { SYLPH_SKETCH; SYLPH_PROFILING } from '../modules/sylph'
include { SOFTWARE_VERSIONS } from '../modules/software_versions'

/* 
 * Main pipeline logic
 */

workflow tofumaapo_sylph {
    main:

        ch_versions = Channel.from([])

        if(!params.step1 && !params.step2){
        // inputs:
            //if(!params.no_qc){
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
            //}
        // SYLPH
            SYLPH_SKETCH(ch_raw_reads)
            ch_versions = ch_versions.mix( SYLPH_SKETCH.out.version.first() )

            if(params.sylph_merge){
                SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches.map{ it ->
                                                                    def metas = "all"
                                                                    return [metas, it[1]]}.groupTuple() )        
                ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version )
            }else{
                SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches )        
                ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version.first() )
            }
        }

        if(params.step1){
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
        }

        if(params.step2){
            if(params.sralist){
                input_sra_list()
                ch_raw_reads = input_sra_list.out.reads
            }else{
                exit 1, "No input in --sralist was declared! This input is needed for step2"
            }

            SYLPH_SKETCH(ch_raw_reads)
            ch_versions = ch_versions.mix( SYLPH_SKETCH.out.version.first() )

            if(params.sylph_merge){
                SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches.map{ it ->
                                                                    def metas = "all"
                                                                    return [metas, it[1]]}.groupTuple() )        
                ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version )
            }else{
                SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches )        
                ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version.first() )
            }
        }

            SOFTWARE_VERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

}