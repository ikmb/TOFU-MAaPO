//
// Check input samplesheet and get reads channel
//

include {
	download_sra
} from '../modules/download_sra.nf'

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
						if(params.assemblymode == "group"){
							def coassemblygroup = row.group //.ifEmpty(exit 1, "Invalid input samplesheet: No group column for coassembly was found")
							if ( coassemblygroup == "null" || coassemblygroup == "") exit 1, "Invalid input samplesheet: No group column for coassembly was found or contains empty fields"
							meta.coassemblygroup = coassemblygroup
						}else if(params.assemblymode == "all"){
							meta.coassemblygroup = 1
						}else if(params.assemblymode == "single"){
							meta.coassemblygroup = meta.id
						}else{ 
							exit 1, "Only allowed modes for coassembly are all, group or single"
						}
						
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
						if(params.assemblymode == "single"){
							meta.coassemblygroup = meta.id
						}else if( params.assemblymode == "all"){
							meta.coassemblygroup = 1
						}else{ 
							exit 1, "Cannot use other modes than single or all for coassembly with this inputmode"
						}
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
						if(params.assemblymode == "group"){
							def coassemblygroup = row.group.ifEmpty(exit 1, "Invalid input samplesheet: No group column for coassembly was found")
							meta.coassemblygroup = coassemblygroup
						}else if(params.assemblymode == "all"){
							meta.coassemblygroup = 1
						}else if(params.assemblymode == "single"){
							meta.coassemblygroup = meta.id
						}else{ 
							exit 1, "Only allowed modes for coassembly are all, group or single"
						}

						if (params.single_end){
							return [meta, [ read1 ] ] 
						}else{
							return [meta, [ read1, read2, read3 ] ]
						}
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
						if(params.assemblymode == "single"){
							meta.coassemblygroup = meta.id
						}else if( params.assemblymode == "all"){
							meta.coassemblygroup = 1
						}else{ 
							exit 1, "Cannot use other modes than single or all for coassembly with this inputmode"
						}

						if (params.single_end){
							return [meta,  row[1] ] 
						}else{  
							return [meta,  row[1] ]
						}
				}
				.set { reads }
		}

	emit:
		reads // channel: [ val(meta), [ 0:read1, 1:read2 ] ] or  [ val(meta), [ 0:readsingle ] ]
}

workflow input_sra {
	main: 
		if(!params.apikey){
			exit 1, "No SRA apikey was declared."
		}else{
			ids = params.sra.split(',').collect()
			Channel.fromSRA(ids, apiKey: params.apikey)
				.collectFile(name: "${params.outdir}/sample_SRA_list.txt", newLine: true)
				.map {row -> 
						def meta = [:]
						meta.id = row[0]

						if(params.assemblymode == "single"){
							meta.coassemblygroup = meta.id
						}else if( params.assemblymode == "all"){
							meta.coassemblygroup = 1
						}else{ 
							exit 1, "Cannot use other modes than single or all for coassembly with this SRA input"
						}

						def singleEnd = row[1].size() == -1
						meta.single_end = singleEnd.toBoolean()
						if( row[1].size() == 3){
							fastq = row[1]
							if (hasExtension(fastq[1], ".fastq.gz") && hasExtension(fastq[2], ".fastq.gz") && hasExtension(fastq[0], ".fastq.gz")){
								//order shall be: meta, forward, reversed, unpaired
								return [meta, [fastq[1], fastq[2], fastq[0]] ]
							}else{
								println "Invalid file paths in triplet ${row[1]}, they do not end with .fastq.gz"
							}
						}else if( row[1].size() == 2){
							fastq = row[1]
							if (hasExtension(fastq[0], ".fastq.gz") && hasExtension(fastq[1], ".fastq.gz")){
								//order shall be: meta, forward, reversed
								return [meta, [fastq[0], fastq[1]] ]
							}else{
								println "Invalid file paths in pair ${row[1]}, they do not end with .fastq.gz"
							}
						}else{
							if (hasExtension(row[1], ".fastq.gz")){
								//order shall be: meta, unpaired
								return [meta,  row[1] ]
							}else{
								println "Invalid file paths in single-end ${row[1]}, it does not end with .fastq.gz"
							}
						}
				}
				.collectFile(name: "${params.outdir}/parsed_sample_list.txt", newLine: true).set { output }
		}
		download_sra(output)
		reads = download_sra.out.reads
	emit:
		reads
}

def hasExtension(it, extension) {
	it.toString().toLowerCase().endsWith(extension.toLowerCase())
}