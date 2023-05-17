//
// Check input samplesheet and get reads channel
//

workflow input_check {
    main:
        if(hasExtension(params.reads, "csv")){
            Channel
                .from(file(params.reads))
                .splitCsv ( header:true, sep:',' )
                .map { row ->
                        def id = row.id
                        def read1 = row.read1 ? file(row.read1, checkIfExists: true) : false
                        if (!read1) exit 1, "Invalid input samplesheet"
                        def read2 = row.read2 ? file(row.read2, checkIfExists: true) : false
                        if (!read2 && !params.single_end) exit 1, "Invalid input samplesheet: Only one read per sample found, but '--single_end' is not set!"
                        def meta = [:]  
                        meta.id = id
                        meta.single_end = params.single_end.toBoolean()
                         if (params.single_end)
                            return [meta, [ read1 ] ] 
                        else  
                            return [meta, [ read1, read2 ] ]
                }
                .set { reads }
        } else {
            Channel
                .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
                .ifEmpty { exit 1, "Cannot find any matching reads in ${params.reads}\nPaths must be in enclosed quotes"}
                .map { row ->
                        def meta = [:]
                        meta.id = row[0]
                        meta.single_end = params.single_end.toBoolean()
                        if (params.single_end)
                            return [meta,  row[1] ] 
                        else  
                            return [meta,  row[1] ]
                }
                .set { reads }
        }

    emit:
        reads // channel: [ val(meta), [ 0:read1, 1:read2 ] ] or  [ val(meta), [ 0:readsingle ] ]
}

workflow input_check_qced {
  
    main:
        if(hasExtension(params.reads, "csv")){
            Channel
                .from(file(params.reads))
                .splitCsv ( header:true, sep:',' )
                .map { row ->
                        def id = row.id

                        def read1 = row.read1 ? file(row.read1, checkIfExists: true) : false
                        if (!read1) exit 1, "Invalid input samplesheet in at least column read1! Is your csv file comma separated?"
                        if (!hasExtension(row.read1, "R1_clean.fastq.gz") ) exit 1, "Invalid read names! Reads need to end with {R1,R2,single}_clean.fastq.gz"

                        def read2 = row.read2 ? file(row.read2, checkIfExists: true) : false
                        if (!read2) exit 1, "Invalid input samplesheet in at least column read2"
                        if ( (!hasExtension(row.read2, "R2_clean.fastq.gz") ) ) exit 1, "Invalid read names! Reads need to end with {R1,R2,single}_clean.fastq.gz"

                        def read3 = row.read3 ? file(row.read3, checkIfExists: true) : false
                        if (!read3) exit 1, "Invalid input samplesheet in at least column read3"
                        if (!hasExtension(row.read3, "single_clean.fastq.gz")) exit 1, "Invalid read names! Reads need to end with {R1,R2,single}_clean.fastq.gz"
                        if (!read3 && !params.single_end) exit 1, "Invalid input samplesheet: Only one or two reads per sample found, either use --single_end for one read per sample or add the fastq file with unpaired reads in a column read3!"

                        def meta = [:]  
                        meta.id = id
                        meta.single_end = params.single_end.toBoolean()
                        if (params.single_end)
                            return [meta, [ read1 ] ] 
                        else  
                            return [meta, [ read1, read2, read3 ] ]
                }
                .set { reads }
        } else {
            Channel
                .fromFilePairs(params.reads, size: params.single_end ? 1 : 3)
                .ifEmpty { exit 1, "Cannot find any matching reads in ${params.reads}\nPaths must be in enclosed quotes"}
                .map { row ->
                        def meta = [:]
                        meta.id = row[0]
                        meta.single_end = params.single_end.toBoolean()
                        if (params.single_end)
                            return [meta,  row[1] ] 
                        else  
                            return [meta,  row[1] ]
                }
                .set { reads }
        }
    
    emit:
        reads // channel: [ val(meta), [ 0:read1, 1:read2 ] ] or  [ val(meta), [ 0:readsingle ] ]
}

workflow input_vamb {
    take: data
    main:
        data.splitCsv ( header:false, sep:',' )
            .map { row ->
                    def meta = [:]  
                    meta.id = row[1]
                    meta.single_end = params.single_end.toBoolean()
                    meta.vamb_group = row[0]
                    if (params.single_end)
                        return [row[1], meta, [ row[3] ] ] 
                    else  
                        return [row[1], meta, [ row[3], row[4], row[5] ] ]
            }
            .set { reads }
    
    emit:
        reads // channel: [ val(meta), [ 0:read1, 1:read2 ] ] or  [ val(meta), [ 0:readsingle ] ]
}


def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}