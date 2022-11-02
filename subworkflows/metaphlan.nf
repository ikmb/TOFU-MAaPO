
/* 
 * Import metaphlan modules 
 */
include {
  PREPARE_METAPHLAN;
  METAPHLAN;
  ABUNDANCE_REL_MERGE;
  ABUNDANCE_ABS_MERGE
  } from '../modules/metaphlan.nf'

/*
 * Metaphlan3 pipeline logic
 */
workflow metaphlan{
    take: data
    main:
        ch_versions = Channel.empty()
        /*
        if (params.metaphlan_db) {
            db_path = file("${params.metaphlan_db_test}")
            if (!db_path.exists()) { exit 1, "Could not find Metaphlan database - please check the path" } 
        }else {exit 1, "No Metaphlan database was specified, aborting..."}
        */
        if(params.updatemetaphlan){
            PREPARE_METAPHLAN()
            ch_readymetaphlan = PREPARE_METAPHLAN.out
        }else{
            ch_readymetaphlan = Channel.of('true')
        }
        METAPHLAN(data, ch_readymetaphlan)
        ch_metaphout = METAPHLAN.out.outputMetaphlan
        ch_versions = ch_versions.mix(METAPHLAN.out.version_metaphlan.first() )
        if(params.metaphlan_analysis_type == "rel_ab_w_read_stats"){
            ABUNDANCE_REL_MERGE(ch_metaphout.collect() )
            ABUNDANCE_ABS_MERGE(ch_metaphout.collect() )
        } else if (params.metaphlan_analysis_type == "rel_ab") {
            ABUNDANCE_REL_MERGE(ch_metaphout.collect() )
        }
    emit:
        versions = ch_versions
}