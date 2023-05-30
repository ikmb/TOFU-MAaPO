/*
 * Import Genome Assembly modules
 */

include { CONTIGS_MAPPING } from '../modules/assembly/bowtie2_contigs_mapping.nf'
include { FILTERCONTIGS } from '../modules/assembly/contig_filter.nf'
include { MEGAHIT } from '../modules/assembly/megahit.nf'
include { checkm } from '../modules/assembly/checkm.nf'
include { MAXBIN2 } from '../modules/assembly/maxbin2.nf'
include { CONCOCT } from '../modules/assembly/concoct.nf'
include { getCountTable } from '../modules/assembly/assembly_util.nf'

include {   GTDBTK; 
            PREPARE_GTDBTK 
} from '../modules/assembly/gtdbtk.nf'

include {   METABAT;
            contigs_to_bins
} from '../modules/assembly/metabat.nf'

include {
    VAMB_CATALOGUE;
    VAMB_CATALOGUE_INDEX;
    VAMB_MAPPING;
    VAMB_COLLECT_DEPTHS;
    VAMB_CONTIGS_SELECTION;
    VAMB;
    group_vamb
    } from '../modules/assembly/vamb.nf'

include { input_vamb } from '../subworkflows/input_check'

include {
    FORMATTING_CONTIG_TO_BIN;
    FORMATTING_CONTIG_TO_BIN_NOVAMB;
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
        ch_versions = Channel.empty()

        if(!params.magscot){
            /*
             * Basic Genome Assembly:
             */
            MEGAHIT(data)
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first() )

            FILTERCONTIGS(MEGAHIT.out.contigs)
            ch_versions = ch_versions.mix(FILTERCONTIGS.out.versions.first() )
            
            ch_filteredcontigs = FILTERCONTIGS.out.contigs

            MINIMAP2_CATALOGUE( ch_filteredcontigs )
            MINIMAP2_CATALOGUE_INDEX( MINIMAP2_CATALOGUE.out.catalogue )


            MINIMAP2_MAPPING( ch_filteredcontigs.join( MINIMAP2_CATALOGUE_INDEX.out.catalogue ) )

            ch_mapping = MINIMAP2_MAPPING.out.maps
            ch_bam = MINIMAP2_MAPPING.out.bam

            METABAT(ch_mapping)
            contigs_to_bins(METABAT.out)
            checkm(METABAT.out)

            if(!params.skip_gtdbtk){
            
                if(params.updategtdbtk){
                    PREPARE_GTDBTK()
                    ch_readygtdbtk = PREPARE_GTDBTK.out.readystate
                }else{
                    ch_readygtdbtk = Channel.of('true')
                }
                GTDBTK(METABAT.out, ch_readygtdbtk)
            }else{
                ch_readygtdbtk = Channel.of('true')
            }
            
            getCountTable(ch_bam)
        
        
            
        }else{
            /*
             * Extended Genome Assembly:
             */
    /*
    * Contigs
    */
            MEGAHIT(data)
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first() )

            FILTERCONTIGS(MEGAHIT.out.contigs)
            ch_versions = ch_versions.mix(FILTERCONTIGS.out.versions.first() )
            
            ch_filteredcontigs = FILTERCONTIGS.out.contigs

            ch_collected_filtered_contigs = FILTERCONTIGS.out.filteredcontig.collect()

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
        if(!params.skip_vamb){
            //create a new csv file to subgroup samples

            ch_allcontigs_table = ch_filteredcontigs.collectFile() { item ->
                [ "contigs_grouped.csv", item[0].id + ',' + item[0].single_end + ',' + item[1] + ',' + item[2] + '\n']
                }

            //new csv file will be read in, we create a file with all fastq files for VAMB_CATALOGUE
            group_vamb(ch_allcontigs_table)


            //get a tuple which sample (and therefor its catalogue) belongs to wich subgroup
            ch_sample_to_vambgroup = group_vamb.out.sample_vambkey
                .splitCsv ( header:false, sep:',' )
                .map { row ->
                        def meta = [:]
                        meta.id = row[0]
                        meta.single_end = row[1].toBoolean()  
                        return [ meta, row[2] ] }//.view()

            //add to tuple meta vamb_group the tuple contigs with reads
            ch_vambgroup_contigs = ch_sample_to_vambgroup.join( ch_filteredcontigs )//.map{row -> tuple(row[1], row[0], row[1], row[2], row[3], row[4])}

            ch_contigs_perkey = group_vamb.out.contigs_perkey
                .splitCsv ( header:false, sep:',' )
                .map { row -> tuple(row[0], row[1]) }

        // Minimap2 Index from all samples
            VAMB_CATALOGUE(ch_contigs_perkey)



            //VAMB_CATALOGUE(ch_collected_filtered_contigs  )
            VAMB_CATALOGUE_INDEX( VAMB_CATALOGUE.out.catalogue )

            //add to ch_vambgroup_contigs the catalogue based on vamb_group
            ch_mapping_vamb_input = ch_vambgroup_contigs.combine( VAMB_CATALOGUE_INDEX.out.catalogue, by: 1 )


            VAMB_MAPPING( ch_mapping_vamb_input )
            
            VAMB_COLLECT_DEPTHS( VAMB_MAPPING.out.counttable.groupTuple()//.collect() 
            )

            VAMB(   VAMB_CATALOGUE_INDEX.out.catalogue_indexfirst.join( VAMB_COLLECT_DEPTHS.out.alldepths )                    
                )

            ch_vambgroup_sampleid = ch_sample_to_vambgroup.map{ row -> tuple(row[1], row[0]) }.combine(VAMB.out.all_samples_clustertable, by: 0)
            

            VAMB_CONTIGS_SELECTION( ch_vambgroup_sampleid )

        }

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

            if(!params.skip_vamb){
                ch_per_sample_contigs_to_bins = VAMB_CONTIGS_SELECTION.out.persample_clustertable.join( contigs_to_bins.out.metabat2_contigs_to_bins ).join( MAXBIN2.out.contigs_to_bin).join( CONCOCT.out.contigs_to_bin )
                FORMATTING_CONTIG_TO_BIN(   ch_per_sample_contigs_to_bins   )
                ch_contig_to_bin = FORMATTING_CONTIG_TO_BIN.out.formatted_contigs_to_bin
            }else{
                ch_per_sample_contigs_to_bins = contigs_to_bins.out.metabat2_contigs_to_bins.join( MAXBIN2.out.contigs_to_bin).join( CONCOCT.out.contigs_to_bin )
                FORMATTING_CONTIG_TO_BIN_NOVAMB(   ch_per_sample_contigs_to_bins   )
                ch_contig_to_bin = FORMATTING_CONTIG_TO_BIN_NOVAMB.out.formatted_contigs_to_bin
            
            }                          
            
            MARKER_IDENT(   ch_mapping.join( ch_contig_to_bin) )
            MAGSCOT( ch_contig_to_bin.join( MARKER_IDENT.out.hmm_output ).join( FILTERCONTIGS.out.magscot_contigs ))
            EXTRACT_REFINED_BINS ( MAGSCOT.out.refined_contigs_to_bins )

    /*
    * Quality Check Workflow
    */
            if(!params.skip_checkm){
                checkm( EXTRACT_REFINED_BINS.out.refined_bins )
            }
            if(!params.skip_gtdbtk){

                if(params.updategtdbtk){
                    PREPARE_GTDBTK()
                    ch_readygtdbtk = PREPARE_GTDBTK.out.readystate
                }else{
                    ch_readygtdbtk = Channel.of('true')
                }
                GTDBTK(EXTRACT_REFINED_BINS.out.refined_bins, ch_readygtdbtk)          

                /*
                * Abundance Table for MAGS
                */



                //BINCOVERAGE_PERSAMPLE( CONTIGS_MAPPING.out.sample_depth.join( MAGSCOT.out.contigs_to_bins_table ).join( GTDBTK.out.taxonomic_table ) )
                BINCOVERAGE_PERSAMPLE( MINIMAP2_MAPPING.out.sample_depth.join( MAGSCOT.out.contigs_to_bins_table ).join( GTDBTK.out.taxonomic_table ) )
            
            }else{
                ch_readygtdbtk = Channel.of('true')
                if(params.updategtdbtk){
                    PREPARE_GTDBTK()
                }
        
            getCountTable(ch_bam)
        
        }
    }
    emit:
        versions = ch_versions
}