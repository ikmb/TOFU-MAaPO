/*
 * Import humann modules
 */
include {
    HUMANN_SE;
    HUMANN_PE;
    JOINgenefamilies;
    JOINpathabundance;
    JOINpathcoverage
    } from '../modules/humann.nf'

/*
 * Humann3 pipeline logic
 */
workflow humann{
    take: data
    main:
        /*
        if (params.metaphlan_db) {
            db_path = file("${params.metaphlan_db_test}")
            if (!db_path.exists()){ exit 1, "Could not find Metaphlan database - please check the path"}
        } else {
	    exit 1, "No Metaphlan database was specified, aborting..."
        }
        //TODO: CHANGE THIS for a path checking solution:
        
        humann_db_path = file("${params.humann_db_test}")
	    if (!humann_db_path.exists() ){ exit 1, "Could not find your HUMAnN DB - please check the path"
        } else if (params.humann) {
        exit 1, "No HUMAnN database was specified, aborting..."
        }
        */
        if(!params.metaphlan){
            PREPARE_METAPHLAN()
        }
        if(params.single_end){
            HUMANN_SE(data)
            JOINgenefamilies(HUMANN_SE.out.genefamilies.collect() )
            JOINpathabundance(HUMANN_SE.out.pathabundance.collect())
            JOINpathcoverage(HUMANN_SE.out.pathcoverage.collect())
        }else{
            HUMANN_PE(data)
            JOINgenefamilies(HUMANN_PE.out.genefamilies.collect() )
            JOINpathabundance(HUMANN_PE.out.pathabundance.collect())
            JOINpathcoverage(HUMANN_PE.out.pathcoverage.collect())
        }
}