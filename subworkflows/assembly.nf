/*
 * Import Genome Assembly modules
 */

//include { CONTIGS_MAPPING } from '../modules/assembly/bowtie2_contigs_mapping.nf' //Currently unused, interchangeable with MINIMAP2_MAPPING
include { FILTERCONTIGS } from '../modules/assembly/contig_filter.nf'
include { MEGAHIT_assembly } from '../modules/assembly/megahit.nf'
include { checkm } from '../modules/assembly/checkm.nf'
include { MAXBIN2 } from '../modules/assembly/maxbin2.nf'
include { CONCOCT } from '../modules/assembly/concoct.nf'
include { SEMIBIN } from '../modules/assembly/semibin.nf'
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

include {
	FORMATTING_CONTIG_TO_BIN;
	MARKER_IDENT;
	MAGSCOT;
	EXTRACT_REFINED_BINS
	} from '../modules/assembly/magscot.nf'

include { 
	BINCOVERAGE_PERSAMPLE;
	MERGE_MAG_ABUNDANCE } from '../modules/assembly/bincoverage_persample.nf'

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

		binner = params.binner ? params.binner.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '').replaceAll('2', '')} : []
		if(params.assembly){log.info "Binner       : ${binner}"}

		ch_contig_bin_list = Channel.empty()

		/*
		* Contigs
		*/
		megahit_coas_input = data.map { it ->
			metas = it[0]
			return[metas.coassemblygroup, it[1]]}.groupTuple(by:0).map{ it -> return[it[0], it[1].flatten()]}.unique()

		MEGAHIT_assembly(megahit_coas_input)
		ch_versions = ch_versions.mix(MEGAHIT_assembly.out.versions.first() )

		filtercontigs_in = MEGAHIT_assembly.out.contigs

		FILTERCONTIGS(filtercontigs_in)
		ch_versions = ch_versions.mix(FILTERCONTIGS.out.versions.first() )

		ch_filteredcontigs = FILTERCONTIGS.out.contigs.combine(data.map { it ->
				metas = it[0]
				return[metas.coassemblygroup, metas, it[1]]}, by:0 ).map { it -> return[it[2], it[1], it[3]]}

		/*
		 * Basic Genome Assembly:
		 */
		MINIMAP2_CATALOGUE( FILTERCONTIGS.out.contigs )
		ch_versions = ch_versions.mix(MINIMAP2_CATALOGUE.out.versions )

		MINIMAP2_CATALOGUE_INDEX( MINIMAP2_CATALOGUE.out.catalogue )
		ch_versions = ch_versions.mix(MINIMAP2_CATALOGUE_INDEX.out.versions )

		ch_minimap2_mapping_input = ch_filteredcontigs.map{it -> meta = it[0]
									return[meta.coassemblygroup, meta, it[1], it[2]]}
									.combine( MINIMAP2_CATALOGUE_INDEX.out.catalogue, by:0 )
									.map{it -> meta = it[1]
									return[meta, it[2], it[3], it[4], it[5]]}

		MINIMAP2_MAPPING( ch_minimap2_mapping_input )

		ch_versions = ch_versions.mix(MINIMAP2_MAPPING.out.versions.first() )

		ch_mapping = MINIMAP2_MAPPING.out.maps
		ch_bam = MINIMAP2_MAPPING.out.bam

		/*
		* METABAT2 Workflow
		*/
		if ( 'metabat' in binner || !params.magscot.toBoolean() ) {
			METABAT(ch_mapping)
			ch_versions = ch_versions.mix(METABAT.out.versions.first() )
		
			contigs_to_bins(METABAT.out.metabatout)
			ch_contig_bin_list = ch_contig_bin_list.mix(contigs_to_bins.out.metabat2_contigs_to_bins)
			if(!params.magscot.toBoolean()){ ch_bins = METABAT.out.metabatout }
		}

		if(params.magscot.toBoolean()){
			/*
			 * Extended Genome Assembly:
			*/

			/*
			* VAMB Workflow
			*/
			if ( 'vamb' in binner ) {

				if(params.assemblymode == "single"){
					//create a new csv file to subgroup samples

					ch_allcontigs_table = ch_filteredcontigs.collectFile() { item ->
						[ "contigs_grouped.csv", item[0].id + ',' + item[0].single_end + ',' +  item[0].coassemblygroup + ',' + item[1] + ',' + item[2] +  '\n']
						}

					//new csv file will be read in, we create a file with all fastq files for VAMB_CATALOGUE
					group_vamb(ch_allcontigs_table)

					//get a tuple which sample (and therefor its catalogue) belongs to wich subgroup
					ch_sample_to_vambgroup = group_vamb.out.sample_vambkey
					.splitCsv ( header:false, sep:',' )
					.map { row ->
							def meta = [:]
							meta.id = row[0]
							meta.coassemblygroup = row[2]
							meta.single_end = row[1].toBoolean()  
							return [ meta, row[3] ] }

					//add to tuple meta vamb_group the tuple contigs with reads
					ch_vambgroup_contigs = ch_sample_to_vambgroup.join( ch_filteredcontigs )//.map{row -> tuple(row[1], row[0], row[1], row[2], row[3], row[4])}

					ch_contigs_perkey = group_vamb.out.contigs_perkey
						.splitCsv ( header:false, sep:',' )
						.map { row -> tuple(row[0], row[1]) }


					vamb_catalogue_in = ch_contigs_perkey


					// Minimap2 Index from all samples
					VAMB_CATALOGUE(vamb_catalogue_in)
					ch_versions = ch_versions.mix(VAMB_CATALOGUE.out.versions.first() )


					//VAMB_CATALOGUE(ch_collected_filtered_contigs  )
					VAMB_CATALOGUE_INDEX( VAMB_CATALOGUE.out.catalogue )
					ch_versions = ch_versions.mix(VAMB_CATALOGUE_INDEX.out.versions.first() )

					//add to ch_vambgroup_contigs the catalogue based on vamb_group
					ch_mapping_vamb_input = ch_vambgroup_contigs.combine( VAMB_CATALOGUE_INDEX.out.catalogue, by: 1 )


					VAMB_MAPPING( ch_mapping_vamb_input )
					ch_versions = ch_versions.mix(VAMB_MAPPING.out.versions.first() )

					VAMB_COLLECT_DEPTHS( VAMB_MAPPING.out.counttable.groupTuple()//.collect() 
					)
					ch_versions = ch_versions.mix(VAMB_COLLECT_DEPTHS.out.versions.first() )

					VAMB(   VAMB_CATALOGUE_INDEX.out.catalogue_indexfirst.join( VAMB_COLLECT_DEPTHS.out.alldepths )                    
						)
					ch_versions = ch_versions.mix(VAMB.out.versions.first() )
					ch_vambgroup_sampleid = ch_sample_to_vambgroup.map{ row -> tuple(row[1], row[0]) }.combine(VAMB.out.all_samples_clustertable, by: 0)
				}else{
					
					VAMB_COLLECT_DEPTHS( MINIMAP2_MAPPING.out.counttable_vamb.groupTuple()//.collect() 
					)
					ch_versions = ch_versions.mix(VAMB_COLLECT_DEPTHS.out.versions.first() )

					VAMB(   MINIMAP2_CATALOGUE_INDEX.out.catalogue_indexfirst.join( VAMB_COLLECT_DEPTHS.out.alldepths )                    
						)
					ch_versions = ch_versions.mix(VAMB.out.versions.first() )


					ch_vambgroup_sampleid = data.map{it -> meta = it[0]
													return[meta.coassemblygroup, meta]}
													.combine(VAMB.out.all_samples_clustertable, by: 0)
				}

				VAMB_CONTIGS_SELECTION( ch_vambgroup_sampleid )
				ch_contig_bin_list = ch_contig_bin_list.mix(VAMB_CONTIGS_SELECTION.out.magscot_contigbinlist)
			}

			/*
			* MAXBIN2 Workflow
			*/
			if ( 'maxbin' in binner ) {
				MAXBIN2( ch_mapping )
				ch_contig_bin_list = ch_contig_bin_list.mix(MAXBIN2.out.magscot_contigbinlist)
				ch_versions = ch_versions.mix(MAXBIN2.out.versions.first() )
			}
			/*
			* CONCOCT Workflow
			*/
			if ( 'concoct' in binner ) {
			CONCOCT( ch_mapping.join( ch_bam ) )
			ch_contig_bin_list = ch_contig_bin_list.mix(CONCOCT.out.magscot_contigbinlist)
			ch_versions = ch_versions.mix(CONCOCT.out.versions.first() )
			}
			/*
			* SemiBin2 Workflow
			*/
			if ( 'semibin' in binner ) {
			SEMIBIN( ch_mapping.join( ch_bam ) )
			ch_contig_bin_list = ch_contig_bin_list.mix(SEMIBIN.out.magscot_contigbinlist)
			ch_versions = ch_versions.mix(SEMIBIN.out.versions.first() )
			}

			/*
			* MAGScoT Workflow
			*/

			ch_per_sample_contigs_to_bins = ch_contig_bin_list.groupTuple()

			FORMATTING_CONTIG_TO_BIN(   ch_per_sample_contigs_to_bins   )
			ch_contig_to_bin = FORMATTING_CONTIG_TO_BIN.out.formatted_contigs_to_bin

			MARKER_IDENT(   ch_mapping.join( ch_contig_to_bin) )
			ch_versions = ch_versions.mix(MARKER_IDENT.out.versions.first() )

			ch_magscot_in = ch_contig_to_bin
									.join( MARKER_IDENT.out.hmm_output )
									.map{it -> meta = it[0]
									return[meta.coassemblygroup, meta, it[1], it[2]]}
									.combine( FILTERCONTIGS.out.magscot_contigs, by:0 )

			MAGSCOT( ch_magscot_in )
			ch_versions = ch_versions.mix(MAGSCOT.out.versions.first() )

			EXTRACT_REFINED_BINS ( MAGSCOT.out.refined_contigs_to_bins )
			ch_versions = ch_versions.mix(EXTRACT_REFINED_BINS.out.versions.first() )

			ch_bins = EXTRACT_REFINED_BINS.out.refined_bins
		}


		/*
		* Quality Check Workflow
		*/
		//TODO: Optimize by only running once for each co-assembly group
		if(!params.skip_checkm){
			checkm( ch_bins )
			ch_versions = ch_versions.mix(checkm.out.versions.first() )
		}
		if(!params.skip_gtdbtk){
			if(params.updategtdbtk){
				PREPARE_GTDBTK()
				ch_readygtdbtk = PREPARE_GTDBTK.out.readystate
			}else{
				ch_readygtdbtk = Channel.of('true')
			}
			GTDBTK(ch_bins, ch_readygtdbtk)
			ch_versions = ch_versions.mix(GTDBTK.out.versions.first() )
			
		}else{
			ch_readygtdbtk = Channel.of('true')
			if(params.updategtdbtk){
				PREPARE_GTDBTK()
			}
		}
		getCountTable(ch_bam)
		ch_versions = ch_versions.mix(getCountTable.out.versions.first() )

		if(params.magscot.toBoolean()){
			/*
			* Abundance Table for MAGS
			*/
			BINCOVERAGE_PERSAMPLE( MINIMAP2_MAPPING.out.sample_depth.join( MAGSCOT.out.contigs_to_bins_table ).join( GTDBTK.out.taxonomic_table ) )
			ch_versions = ch_versions.mix(BINCOVERAGE_PERSAMPLE.out.versions.first() )

			MERGE_MAG_ABUNDANCE(BINCOVERAGE_PERSAMPLE.out.abundancetable.collect() )
			ch_versions = ch_versions.mix(MERGE_MAG_ABUNDANCE.out.versions.first() )

		}


	emit:
		versions = ch_versions
}