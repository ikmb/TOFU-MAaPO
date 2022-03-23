/* 
 * Import kraken modules 
 */
include {
  KRAKEN2;
//  KRAKEN2_SE;
//  KRAKEN2_PE;
  KRAKEN2MPA;
  KRAKEN2YAML;
  KRAKENMERGEREPORTS;
  KRAKENMPAMERGE;
  BRACKEN;
  BRACKENMERGE
  } from '../modules/kraken.nf' 
  
/*
 * Kraken2 pipeline logic
 */
workflow kraken{
    take: data
    main:
    /*
        if(params.single_end){
            KRAKEN2_SE(data)
            KRAKEN2MPA(KRAKEN2_SE.out.krakenreport)
            KRAKEN2YAML(KRAKEN2_SE.out.krakenreport.collect()  )
            KRAKENMERGEREPORTS(KRAKEN2_SE.out.krakenreport.collect()   )
            KRAKENMPAMERGE(KRAKEN2MPA.out.krakenmpa.collect()  )
                if(params.bracken){
                    BRACKEN(KRAKEN2_SE.out.brackeninput)
                    BRACKENMERGE(BRACKEN.out.collect()   )
                }
        }else{
            KRAKEN2_PE(data)
            KRAKEN2MPA(KRAKEN2_PE.out.krakenreport)
            KRAKEN2YAML(KRAKEN2_PE.out.krakenreport.collect()  )
            KRAKENMERGEREPORTS(KRAKEN2_PE.out.krakenreport.collect()   )
            KRAKENMPAMERGE(KRAKEN2MPA.out.krakenmpa.collect()  )
                if(params.bracken){
                    BRACKEN(KRAKEN2_PE.out.brackeninput)
                    BRACKENMERGE(BRACKEN.out.collect()   )
                }
        }
    */
    KRAKEN2(data)
            KRAKEN2MPA(KRAKEN2.out.krakenreport)
            KRAKEN2YAML(KRAKEN2.out.krakenreport.collect()  )
            KRAKENMERGEREPORTS(KRAKEN2.out.krakenreport.collect()   )
            KRAKENMPAMERGE(KRAKEN2MPA.out.krakenmpa.collect()  )
                if(params.bracken){
                    BRACKEN(KRAKEN2.out.brackeninput)
                    BRACKENMERGE(BRACKEN.out.collect()   )
                }
    emit:
        kraken_data = KRAKEN2YAML.out
}