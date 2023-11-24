include { SALMON; SALMON_merge } from '../modules/salmon'

workflow salmon{
	take: data
	main:
		ch_versions = Channel.empty()

        SALMON(data)
        SALMON_merge(SALMON.out.salmonout.collect() )

        ch_versions = ch_versions.mix( SALMON.out.version.first() )
		ch_versions = ch_versions.mix( SALMON_merge.out.version )

	emit:
		salmon_data = SALMON_merge.out.salmonmerge
		versions = ch_versions
}