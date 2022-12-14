/*
 * Import Genome Assembly modules
 */

include { CONTIGS_MAPPING } from '../modules/assembly/bowtie2_contigs_mapping.nf'
include { filtercontigs } from '../modules/assembly/contig_filter.nf'
include { MEGAHIT } from '../modules/assembly/megahit.nf'
include { GTDBTK } from '../modules/assembly/gtdbtk.nf'
include { checkm } from '../modules/assembly/checkm.nf'
include { MAXBIN2 } from '../modules/assembly/maxbin2.nf'
include { CONCOCT } from '../modules/assembly/concoct.nf'
include { getCountTable } from '../modules/assembly/assembly_util.nf'

include { METABAT;
          contigs_to_bins
        } from '../modules/assembly/metabat.nf'

include {
    VAMB_CATALOGUE;
    VAMB_CATALOGUE_INDEX;
    VAMB_MAPPING;
    VAMB_COLLECT_DEPTHS;
    VAMB_CONTIGS_SELECTION;
    VAMB
    } from '../modules/assembly/vamb.nf'

include {
    FORMATTING_CONTIG_TO_BIN;
    MARKER_IDENT;
    MAGSCOT;
    EXTRACT_REFINED_BINS
    } from '../modules/assembly/magscot.nf'

include { BINCOVERAGE_PERSAMPLE } from '../modules/assembly/bincoverage_persample.nf'

include { 
        MINIMAP2_CATALOGUE;
        MINIMAP2_CATALOGUE_INDEX;
        MINIMAP2_MAPPING } from '../modules/assembly/minimap2_mapping.nf'

/*
 * Genome Assembly pipeline logic
 */
workflow assembly{
    take: data
    main:
        if(!params.magscot){
            MEGAHIT(data)
            filtercontigs(MEGAHIT.out.contigs)
            
            ch_filteredcontigs = filtercontigs.out.contigs

            //CONTIGS_MAPPING(ch_filteredcontigs)


            MINIMAP2_CATALOGUE(ch_collected_filtered_contigs  )
            MINIMAP2_CATALOGUE_INDEX( MINIMAP2_CATALOGUE.out.catalogue )


            MINIMAP2_MAPPING(   ch_filteredcontigs, 
                                MINIMAP2_CATALOGUE_INDEX.out.catalogue,
                              )

            ch_mapping = MINIMAP2_MAPPING.out.maps
            ch_bam = MINIMAP2_MAPPING.out.bam

            METABAT(ch_mapping)
            contigs_to_bins(METABAT.out)
            checkm(METABAT.out)

            if(!params.skip_gtdbtk){
                GTDBTK(METABAT.out)
            }
        
            getCountTable(ch_bam)
        }else{

    /*
    * Contigs
    */
            MEGAHIT(data)
            filtercontigs(MEGAHIT.out.contigs)
            
            ch_filteredcontigs = filtercontigs.out.contigs

            //TODO: Create batches of e.g. 100 samples for vamb; needs same samples in catalogue creation and VAMB_COLLECT_DEPTHS
            //ch_collected_filtered_contigs = filtercontigs.out.filteredcontig.collect().flatten().buffer( size: 100, remainder: true )

            ch_collected_filtered_contigs = filtercontigs.out.filteredcontig.collect()

        // Minimap2 Index from each individual sample
            MINIMAP2_CATALOGUE(ch_filteredcontigs  )
            MINIMAP2_CATALOGUE_INDEX( MINIMAP2_CATALOGUE.out.catalogue )

            MINIMAP2_MAPPING( ch_filteredcontigs.join( MINIMAP2_CATALOGUE_INDEX.out.catalogue ) )
                              
            ch_mapping = MINIMAP2_MAPPING.out.maps
            ch_bam = MINIMAP2_MAPPING.out.bam

    /*
    * METABAT2 Workflow
    */
            METABAT(ch_mapping)
            contigs_to_bins(METABAT.out)

    /*
    * VAMB Workflow
    */
//TODO: MAKE VAMB OPTIONAL
        // Minimap2 Index from all samples

            VAMB_CATALOGUE(ch_collected_filtered_contigs  )
            VAMB_CATALOGUE_INDEX( VAMB_CATALOGUE.out.catalogue )


            VAMB_MAPPING(   ch_filteredcontigs, 
                            VAMB_CATALOGUE_INDEX.out.catalogue
                         )
            
            VAMB_COLLECT_DEPTHS( VAMB_MAPPING.out.counttable.collect() )

            VAMB(   VAMB_CATALOGUE_INDEX.out.catalogue,
                    VAMB_COLLECT_DEPTHS.out.alldepths
                 )
            
            VAMB_CONTIGS_SELECTION( VAMB.out.all_samples_clustertable,
                                    MINIMAP2_MAPPING.out.sampleid
                                    )
            
    /*
    * MAXBIN2 Workflow
    */
            MAXBIN2( ch_mapping )

    /*
    * CONCOCT Workflow
    */
            CONCOCT( ch_mapping.join( ch_bam ) )


    /*
    * MAGScoT Workflow
    */

            ch_per_sample_contigs_to_bins = VAMB_CONTIGS_SELECTION.out.persample_clustertable.join( contigs_to_bins.out.metabat2_contigs_to_bins ).join( MAXBIN2.out.contigs_to_bin).join( CONCOCT.out.contigs_to_bin )
                          
            FORMATTING_CONTIG_TO_BIN(   ch_per_sample_contigs_to_bins   )
            MARKER_IDENT(   ch_mapping.join( FORMATTING_CONTIG_TO_BIN.out.formatted_contigs_to_bin))
            MAGSCOT( FORMATTING_CONTIG_TO_BIN.out.formatted_contigs_to_bin.join( MARKER_IDENT.out.hmm_output ).join( filtercontigs.out.magscot_contigs ))
            EXTRACT_REFINED_BINS ( MAGSCOT.out.refined_contigs_to_bins )
    /*
    * Quality Check Workflow
    */

            checkm( EXTRACT_REFINED_BINS.out.refined_bins )

            if(!params.skip_gtdbtk){
                GTDBTK( EXTRACT_REFINED_BINS.out.refined_bins_folder )
                /*
                * Abundance Table for MAGS
                */
                //BINCOVERAGE_PERSAMPLE( CONTIGS_MAPPING.out.sample_depth.join( MAGSCOT.out.contigs_to_bins_table ).join( GTDBTK.out.taxonomic_table ) )
                BINCOVERAGE_PERSAMPLE( MINIMAP2_MAPPING.out.sample_depth.join( MAGSCOT.out.contigs_to_bins_table ).join( GTDBTK.out.taxonomic_table ) )

            
            }
        
            getCountTable(ch_bam)

        }

        
        
}