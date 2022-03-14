/*
 * Import Genome Assembly modules
 */
include {
    MEGAHIT_PE;
    MEGAHIT_SE;
    MAPPING_PE;
    MAPPING_SE;
    METABAT;
    filtercontigs_SE;
    filtercontigs_PE;
    contigs_to_bins;
    checkm_all_bins;
    GTDBTK;
    getCountTable
    } from '../modules/assembly.nf'

/*
 * Genome Assembly pipeline logic
 */
workflow assembly{
    take: data
    main:
        if(!params.single_end){
            MEGAHIT_PE(data)
            filtercontigs_PE(MEGAHIT_PE.out.contigs)
            
            ch_filteredcontigs = filtercontigs_PE.out.contigs
            MAPPING_PE(ch_filteredcontigs)
            ch_mapping = MAPPING_PE.out.maps
            ch_counttable = MAPPING_PE.out.counttable
        }else{
            MEGAHIT_SE(data)
            filtercontigs_SE(MEGAHIT_SE.out.contigs)

            ch_filteredcontigs = filtercontigs_SE.out.contigs
            MAPPING_SE(ch_filteredcontigs)
            ch_mapping = MAPPING_SE.out.maps
            ch_counttable = MAPPING_SE.out.counttable
        }
        
        METABAT(ch_mapping)
        contigs_to_bins(METABAT.out)
        checkm_all_bins(METABAT.out)
        GTDBTK(METABAT.out)
        getCountTable(ch_counttable)
}

