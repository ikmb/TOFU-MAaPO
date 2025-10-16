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

/*
Workflow:
1. Assembly with megahit
2. Filter contigs by length (2000bp default)
3. Mapping of a catalogue and indexing with minimap2
4. Mapping of reads to the catalgoue, coverage and depth tables + bam files
5. Binning with multiple tools: metabat, vamb, maxbin, concoct, semibin
6. Refinement of multiple binning tool outputs with magscot to a single bin set per sample
7. Quality check with CheckM and taxonomical assignment with GTDB-Tk

*/
workflow assembly{
	take: data //tuple of (meta,reads)
	main:
		//collecting softerware versions in a single channel
		ch_versions = Channel.empty()

		//parsing the chosen set of binning tools
		binner = params.binner ? params.binner.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '').replaceAll('2', '')} : []
		if(params.assembly){log.info "Binner       : ${binner}"}

		//channel for all retained bins
		ch_contig_bin_list = Channel.empty()

		/*
		* Contigs
		*/
		if(params.assemblymode == "single"){ //assemble each sample separately
			megahit_coas_input = data.unique()
				.map { it ->
					metas = it[0]
					return[metas.coassemblygroup, it[1].flatten()]}
		}else{ //assemble per coassemblygroup, collect all reads from samples belonging to the same group and create tuple for assembly
			megahit_coas_input = data.unique()
				.map { it ->
					metas = it[0]
					return[metas.coassemblygroup, it[1]]}
				.flatMap { id, reads -> 
					reads.collect { read -> [id, read] }
				}
				.groupTuple()
		}
		//run megahit assembly per coassemblygroup
		MEGAHIT_assembly(megahit_coas_input)
		ch_versions = ch_versions.mix(MEGAHIT_assembly.out.versions.first() )

		//filter contigs by length
		filtercontigs_in = MEGAHIT_assembly.out.contigs

		FILTERCONTIGS(filtercontigs_in)
		ch_versions = ch_versions.mix(FILTERCONTIGS.out.versions.first() )

		//align per-sample metadata to filtered contigs 
		ch_filteredcontigs = FILTERCONTIGS.out.contigs.combine(data.map { it ->
				metas = it[0]
				return[metas.coassemblygroup, metas, it[1]]}, by:0 ).map { it -> return[it[2], it[1], it[3]]}

		/*
		 * Per sample binning: Mapping and indexing of a catalogue per coassemblygroup and map reads to it
		 */
		MINIMAP2_CATALOGUE( FILTERCONTIGS.out.contigs )
		ch_versions = ch_versions.mix(MINIMAP2_CATALOGUE.out.versions.first() )

		MINIMAP2_CATALOGUE_INDEX( MINIMAP2_CATALOGUE.out.catalogue )
		ch_versions = ch_versions.mix(MINIMAP2_CATALOGUE_INDEX.out.versions.first() )

		ch_minimap2_mapping_input = ch_filteredcontigs.map{it -> meta = it[0]
									return[meta.coassemblygroup, meta, it[1], it[2]]}
									.combine( MINIMAP2_CATALOGUE_INDEX.out.catalogue, by:0 )
									.map{it -> meta = it[1]
									return[meta, it[2], it[3], it[4], it[5]]}

		MINIMAP2_MAPPING( ch_minimap2_mapping_input )

		ch_versions = ch_versions.mix(MINIMAP2_MAPPING.out.versions.first() )
		
		//outputs used for all selected per sample binning tools:
		ch_mapping = MINIMAP2_MAPPING.out.maps
		ch_bam = MINIMAP2_MAPPING.out.bam

		/*
		* METABAT2 Workflow
		*/
		//If magscot is disabled, MetaBAT bins are the final bin output.
		if ( 'metabat' in binner || !params.magscot.toBoolean() ) {
			METABAT(ch_mapping)
			ch_versions = ch_versions.mix(METABAT.out.versions.first() )
		
			contigs_to_bins(METABAT.out.metabatout)
			ch_contig_bin_list = ch_contig_bin_list.mix(contigs_to_bins.out.metabat2_contigs_to_bins)
			//if magscot is disabled, directly use metabat bins for final quality check and taxnomic assignment
			if(!params.magscot.toBoolean()){ ch_bins = METABAT.out.metabatout }
		}

		if(params.magscot.toBoolean()){
			/*
			 * Extended Genome Assembly:
			 * MAGScoT bin refinement allows for combining multiple binning tools outputs into a unified set of high-quality bins
			*/

			/*
			* VAMB Workflow
			*/
			//Co-binning with VAMB (we will require to group multiple samples together for good results, user adjustable group sizes (100 as default))
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


					// Minimap2 catalogue for every co-binning group
					VAMB_CATALOGUE(vamb_catalogue_in)
					ch_versions = ch_versions.mix(VAMB_CATALOGUE.out.versions.first() )


					//VAMB_CATALOGUE(ch_collected_filtered_contigs  )
					VAMB_CATALOGUE_INDEX( VAMB_CATALOGUE.out.catalogue )
					ch_versions = ch_versions.mix(VAMB_CATALOGUE_INDEX.out.versions.first() )

					//add to ch_vambgroup_contigs the catalogue based on vamb_group
					ch_mapping_vamb_input = ch_vambgroup_contigs.combine( VAMB_CATALOGUE_INDEX.out.catalogue, by: 1 )

					//Map reads to the co-binning catalogue and retain depths
					VAMB_MAPPING( ch_mapping_vamb_input )
					ch_versions = ch_versions.mix(VAMB_MAPPING.out.versions.first() )

					VAMB_COLLECT_DEPTHS( VAMB_MAPPING.out.counttable.groupTuple()//.collect() 
					)
					ch_versions = ch_versions.mix(VAMB_COLLECT_DEPTHS.out.versions.first() )
					//Main VAMB bin clustering process
					VAMB(   VAMB_CATALOGUE_INDEX.out.catalogue_indexfirst.join( VAMB_COLLECT_DEPTHS.out.alldepths )                    
						)
					ch_versions = ch_versions.mix(VAMB.out.versions.first() )

					//map vamb clusters back to orignal samples
					ch_vambgroup_sampleid = ch_sample_to_vambgroup.map{ row -> tuple(row[1], row[0]) }.combine(VAMB.out.all_samples_clustertable, by: 0)
				}else{
					//when co-assembly was performed, use earlier produced minimap2 counttables for depths directly
					VAMB_COLLECT_DEPTHS( MINIMAP2_MAPPING.out.counttable_vamb.groupTuple()//.collect() 
					)
					ch_versions = ch_versions.mix(VAMB_COLLECT_DEPTHS.out.versions.first() )
					//Main VAMB bin clustering process
					VAMB(   MINIMAP2_CATALOGUE_INDEX.out.catalogue_indexfirst.join( VAMB_COLLECT_DEPTHS.out.alldepths )                    
						)
					ch_versions = ch_versions.mix(VAMB.out.versions.first() )

					//map vamb clusters back to orignal samples
					ch_vambgroup_sampleid = data.map{it -> meta = it[0]
													return[meta.coassemblygroup, meta]}
													.combine(VAMB.out.all_samples_clustertable, by: 0)
				}
				//convert vamb clusters to a contig-bin assignment table
				VAMB_CONTIGS_SELECTION( ch_vambgroup_sampleid )
				ch_contig_bin_list = ch_contig_bin_list.mix(VAMB_CONTIGS_SELECTION.out.magscot_contigbinlist)
			}

			/*
			* MAXBIN2 Workflow
			*/
			if ( 'maxbin' in binner ) {
				MAXBIN2( ch_mapping ) //output is a contig-bin assingment table
				ch_contig_bin_list = ch_contig_bin_list.mix(MAXBIN2.out.magscot_contigbinlist)
				ch_versions = ch_versions.mix(MAXBIN2.out.versions.first() )
			}
			/*
			* CONCOCT Workflow
			*/
			//concoct requires BAM files as input
			if ( 'concoct' in binner ) {
			CONCOCT( ch_mapping.join( ch_bam ) ) //output is a contig-bin assingment table
			ch_contig_bin_list = ch_contig_bin_list.mix(CONCOCT.out.magscot_contigbinlist)
			ch_versions = ch_versions.mix(CONCOCT.out.versions.first() )
			}
			/*
			* SemiBin2 Workflow
			*/
			if ( 'semibin' in binner ) {
			//semibin requires both BAM files and coverage as input
			SEMIBIN( ch_mapping.join( ch_bam ) ) //output is a contig-bin assingment table
			ch_contig_bin_list = ch_contig_bin_list.mix(SEMIBIN.out.magscot_contigbinlist)
			ch_versions = ch_versions.mix(SEMIBIN.out.versions.first() )
			}

			/*
			* MAGScoT Workflow
			*/
			//Aggregate all contig-bin assingment tables for each sample
			ch_per_sample_contigs_to_bins = ch_contig_bin_list.groupTuple()

			FORMATTING_CONTIG_TO_BIN(   ch_per_sample_contigs_to_bins   )
			ch_contig_to_bin = FORMATTING_CONTIG_TO_BIN.out.formatted_contigs_to_bin
			//detect marker genes
			MARKER_IDENT(   ch_mapping.join( ch_contig_to_bin) )
			ch_versions = ch_versions.mix(MARKER_IDENT.out.versions.first() )

			//magscot input is contig-bin assignment table + markers + filtered contigs
			ch_magscot_in = ch_contig_to_bin
									.join( MARKER_IDENT.out.hmm_output )
									.map{it -> meta = it[0]
									return[meta.coassemblygroup, meta, it[1], it[2]]}
									.combine( FILTERCONTIGS.out.magscot_contigs, by:0 )
			//Main MAGScoT bin refinement process, output is a contig-bin assignment table
			MAGSCOT( ch_magscot_in )
			ch_versions = ch_versions.mix(MAGSCOT.out.versions.first() )
			//create final bins
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
			ch_versions = ch_versions.mix(MERGE_MAG_ABUNDANCE.out.versions )

		}


	emit:
		versions = ch_versions
}