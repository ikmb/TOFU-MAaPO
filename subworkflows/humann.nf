/*
 * Import humann modules
 */
include {
    HUMANN;
    PREPARE_HUMANN;
    JOINgenefamilies;
    JOINpathabundance;
    JOINpathcoverage
    } from '../modules/humann.nf'

include {
    PREPARE_METAPHLAN
} from '../modules/metaphlan.nf'
/*
 * Humann3 pipeline logic
 */
workflow humann{
    take: data
    main:

        if(!params.metaphlan){
            if(params.updatemetaphlan){
                PREPARE_METAPHLAN()
                ch_readymetaphlan = PREPARE_METAPHLAN.out.readystate
            }else{
                ch_readymetaphlan = Channel.of('true')
            }
        }else{
                ch_readymetaphlan = Channel.of('true')
        }

        if(params.updatehumann){
            PREPARE_HUMANN()
            ch_readyhumann = PREPARE_HUMANN.out.readystate
        }else{
                ch_readyhumann = Channel.of('true')
        }

        HUMANN(data,ch_readymetaphlan, ch_readyhumann)
        JOINgenefamilies(HUMANN.out.genefamilies.collect() )
        JOINpathabundance(HUMANN.out.pathabundance.collect())
        JOINpathcoverage(HUMANN.out.pathcoverage.collect())

}