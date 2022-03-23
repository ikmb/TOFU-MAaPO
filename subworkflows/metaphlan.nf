
/* 
 * Import metaphlan modules 
 */
include {
  PREPARE_METAPHLAN;
  METAPHLAN;
  //METAPHLAN_SE;
  //METAPHLAN_PE;
  ABUNDANCE_REL_MERGE;
  ABUNDANCE_ABS_MERGE
  } from '../modules/metaphlan.nf'

/*
 * Metaphlan3 pipeline logic
 */
workflow metaphlan{
    take: data
    main:
        /*
        if (params.metaphlan_db) {
            db_path = file("${params.metaphlan_db_test}")
            if (!db_path.exists()) { exit 1, "Could not find Metaphlan database - please check the path" } 
        }else {exit 1, "No Metaphlan database was specified, aborting..."}
        */
        if(params.updatemetaphlan){
            PREPARE_METAPHLAN()
        }
        /*
        if(params.single_end){
            METAPHLAN_SE(data)
            ch_metaphout = METAPHLAN_SE.out.outputMetaphlan
        }else{
            METAPHLAN_PE(data)
            ch_metaphout = METAPHLAN_PE.out.outputMetaphlan
        }
        */
        METAPHLAN(data)
        ch_metaphout = METAPHLAN.out.outputMetaphlan

        if(params.metaphlan_analysis_type == "rel_ab_w_read_stats"){
            ABUNDANCE_REL_MERGE(ch_metaphout.collect() )
            ABUNDANCE_ABS_MERGE(ch_metaphout.collect() )
        } else if (params.metaphlan_analysis_type == "rel_ab") {
            ABUNDANCE_REL_MERGE(ch_metaphout.collect() )
        }
}