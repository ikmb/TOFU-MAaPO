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
        ch_versions = Channel.empty()

        KRAKEN2(data)
        KRAKEN2MPA(KRAKEN2.out.krakenreport)
        KRAKEN2YAML(KRAKEN2.out.krakenreport.collect()  )
        KRAKENMERGEREPORTS(KRAKEN2.out.krakenreport.collect()   )
        KRAKENMPAMERGE(KRAKEN2MPA.out.krakenmpa.collect()  )

        ch_versions = ch_versions.mix( KRAKEN2.out.version.first() )
        ch_versions = ch_versions.mix( KRAKEN2MPA.out.version.first() )

        ch_versions = ch_versions.mix( KRAKEN2YAML.out.version )
        ch_versions = ch_versions.mix( KRAKENMERGEREPORTS.out.version )
        ch_versions = ch_versions.mix( KRAKENMPAMERGE.out.version )


        if(params.bracken){
            BRACKEN(KRAKEN2.out.brackeninput)
            BRACKENMERGE(BRACKEN.out.brackenoutput.collect()   )
            ch_versions = ch_versions.mix(BRACKEN.out.version.first() )
            ch_versions = ch_versions.mix(BRACKENMERGE.out.version )
        }

    emit:

        kraken_data = KRAKEN2YAML.out.kraken2yaml
        versions = ch_versions
}