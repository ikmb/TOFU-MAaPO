include { SYLPH_SKETCH; SYLPH_PROFILING } from '../modules/sylph'

workflow sylph{
	take: data
	main:
		ch_versions = Channel.empty()
		
        SYLPH_SKETCH(data)
		ch_versions = ch_versions.mix( SYLPH_SKETCH.out.version.first() )


		if(params.sylph_merge){
			SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches.map{ it ->
																def metas = "all"
																return [metas, it[1]]}.groupTuple().combine(sylph_database) )
			ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version )
		}else{
			SYLPH_PROFILING(SYLPH_SKETCH.out.sylph_sketches.combine(sylph_database) )      
			ch_versions = ch_versions.mix( SYLPH_PROFILING.out.version.first() )
		}
	emit:
		salmon_data = SYLPH_PROFILING.out.results.collect()
		versions = ch_versions
}