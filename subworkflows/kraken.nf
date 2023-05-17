/* 
 * Import kraken modules 
 */
include {
  KRAKEN2;
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