/*
 * Import humann modules
 */
include {
    HUMANN;
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

        if(!params.metaphlan){
            PREPARE_METAPHLAN()
        }

        HUMANN(data)
        JOINgenefamilies(HUMANN.out.genefamilies.collect() )
        JOINpathabundance(HUMANN.out.pathabundance.collect())
        JOINpathcoverage(HUMANN.out.pathcoverage.collect())

}