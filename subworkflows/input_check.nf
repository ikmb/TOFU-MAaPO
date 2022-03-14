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
        reads // channel: [ val(meta), [ read1, read2 ] ] or  [ val(meta), [ read1 ] ]
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}