process MIN_READS {
    tag "${meta.id}"
    label 'seqkit'
    label 'short_run'

    publishDir "${params.outdir}",
        mode: 'copy',
        pattern: "*_failed.txt",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(meta), path(reads, name: 'input/*')
    val minimum
    val stage

    output:
    tuple val(meta), path("*.fastq.gz"), optional: true, emit: fastq
    tuple val(meta), path(sample_failed), optional: true, emit: failed_log

    script:
    def files = reads instanceof List ? reads : [reads]

    def input_files = files.join(' ')

    def links = files
        .collect { file ->
            def basename = file.fileName.name
            "ln -s ${file} ${basename}"
        }
        .join('\n')

    def sample_failed = "${meta.id}_failed.txt"

    """
    count=\$(seqkit stats -T ${input_files} |
        awk 'NR > 1 { n += \$4 } END { print n+0 }')

    echo 'Sample "${meta.id}" in stage ${stage} contains '\${count}' reads' >&2

    if (( count >= ${minimum} )); then
        ${links}
    else
        echo "Skipping ${meta.id}: ${stage} reads \$count < ${minimum}" >&2
        echo "Skipping ${meta.id}: ${stage} reads \$count < ${minimum}" > ${sample_failed}
    fi
    """
}