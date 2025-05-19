//
// Check input samplesheet and get reads channel
//

include {
	download_sra;
	download_files
} from '../modules/download.nf'

include {
	query_metadata
} from '../modules/entrez.nf'

include { COLLECTOR } from '../modules/QC/collect.nf'

workflow input_check {
	main:
		if(hasExtension(params.reads, "csv")){
			Channel
				.from(file(params.reads))
				.splitCsv ( header:true, sep:',' )
				.map { row ->
						
						def id = row.id ? row.id : false
						if(!id){ exit 1, "Empty id in csv-input found" }
						def meta = [:]  
						meta.id = id
						performdownload = false

						def read1 = row.read1 ? row.read1 : false
						def read2 = row.read2 ? row.read2 : false
						def read3 = row.read3 ? row.read3 : false
						def readsize = [ read1, read2, read3 ].count{ it }

						meta.single_end = readsize == 1 ? true : false

						if(!hasinternetprefix(read1) ){
							// FILE IS LOCALLY AVAILABLE
							if (!FileCheck.checkoutfile("$read1")) { exit 1, "Within the input csv the file $read1 in column read1 does not exist for $id"}
							if (readsize >= 2){ if (!FileCheck.checkoutfile("$read2")) { exit 1, "Within the input csv the file $read2 in column read2 does not exist for $id"} }
							if (readsize == 3){ if (!FileCheck.checkoutfile("$read3")) { exit 1, "Within the input csv the file $read3 in column read3 does not exist for $id"} }
							if (!read1) exit 1, "Invalid input samplesheet in at least column read1! Is your csv file comma separated?"
							if (!(hasExtension(read1, ".fastq.gz") || hasExtension(read1, ".fq.gz")) ) exit 1, "Invalid file $read1 ! Reads need to end with .fastq.gz"

							if ( readsize >= 2 ){
								if (!read2) exit 1, "Invalid input samplesheet in at least column read2"
								if (!FileCheck.checkoutfile("$read2")) { exit 1, "Within the input csv the file $read2 in column read2 does not exist for $id"}
								if (!(hasExtension(row.read2, ".fastq.gz") || hasExtension(row.read2, ".fq.gz")) ) exit 1, "Invalid file $read2 ! Reads need to end with .fastq.gz"
							}
							if (readsize == 3 ){
								if (!read3) exit 1, "Invalid input samplesheet in at least column read3"
								if (!FileCheck.checkoutfile("$read3")) { exit 1, "Within the input csv the file $read3 in column read3 does not exist for $id"}
								if (!(hasExtension(row.read3, ".fastq.gz") || hasExtension(row.read3, ".fq.gz")) ) exit 1, "Invalid file $read3 ! Reads need to end with .fastq.gz"
							}
						}else{
							// FILE IS ONLINE AVAILABLE
							performdownload = true
							//TODO: Ceck if file is valid and available
							if (!(hasExtension(read1, ".fastq.gz") || hasExtension(read1, ".fq.gz")) ) exit 1, "Invalid file $read1 ! Reads need to end with .fastq.gz"
							if ( readsize >= 2 ){
								if (!read2) exit 1, "Invalid input samplesheet in at least column read2"
								if (!(hasExtension(row.read2, ".fastq.gz") || hasExtension(row.read2, ".fq.gz")) ) exit 1, "Invalid file $read2 ! Reads need to end with .fastq.gz"
							}
							if (readsize == 3 ){
								if (!read3) exit 1, "Invalid input samplesheet in at least column read3"
								if (!(hasExtension(row.read3, ".fastq.gz") || hasExtension(row.read3, ".fq.gz")) ) exit 1, "Invalid file $read3 ! Reads need to end with .fastq.gz"
							}
						}

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
						
						if (readsize == 1){
								return [meta, performdownload, [ read1 ] ] 
						}else if (readsize == 2){
								return [meta, performdownload, [ read1, read2 ] ] 
						}else if (readsize == 3){
								return [meta, performdownload, [ read1, read2, read3 ] ] 
						}
				}.branch { it ->
						download: it[1] == true
						local: it[1] == false
					}.set { rawoutput }

				download_files(rawoutput.download.map{it -> return [it[0], it[2]]})
				ch_mix = Channel.empty().mix(rawoutput.local.map{it -> return [it[0], it[2]]}).mix(download_files.out.reads)
				inreads = ch_mix
		} else {
			Channel
				.fromFilePairs(params.reads, size: params.single_end ? 1 : ( params.no_qc ? (params.only_paired ? 2 : 3) : 2 ) ) 
				.ifEmpty { exit 1, 
					( !params.no_qc ? "Cannot find any matching reads in ${params.reads}\nPaths must be in enclosed quotes.\nIf you have single-end reads, make sure to add the parameter --single_end." :
					"Cannot find three quality controlled matching reads (R1, R2 and unpaired) reads in ${params.reads}!\nIf you want to include only paired reads, use the parameter --only_paired.\nIf you have single-end reads, make sure to add the parameter --single_end.")}
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
				.set { inreads }
		}

        if(params.no_qc){
            COLLECTOR( inreads )
		    reads = COLLECTOR.out 
        }else{
            reads = inreads
        }

	emit:
		reads // channel: [ val(meta), [ 0:read1, 1:read2 ] ] or  [ val(meta), [ 0:readsingle ] ]
}

workflow input_sra {
	main: 
		if(!params.apikey){
			exit 1, "No SRA apikey was declared."
		}else{
			ids = params.sra.split(',').collect{it.trim()}
			Channel.fromSRA(ids, apiKey: params.apikey)
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
				.set { raw_rawoutput }
		}
		if(!params.exact_matches){
            rawoutput = raw_rawoutput
        }else{
            rawoutput = raw_rawoutput.filter { meta, files ->
                def sampleid = meta.id
                sampleid in ids
            }
        }
		ch_collectedinput = rawoutput.collectFile(storeDir: "${params.outdir}", name: "parsed_sample_list.csv" ) { item ->
						item[0].id + ',' + item[0].single_end + ',' +  item[0].coassemblygroup + ',' + item[1] + '\n'
						}
		if(!params.step1){

			ch_collectedinput
			.splitCsv ( header:false, sep:',' )
			.map {row -> 
					def read1 = row[3] ? row[3].replaceAll(/\[/, "").replaceAll(/\]/, "").replaceAll(/ /, "") : null
					if (!read1) exit 1, "Invalid input samplesheet in at least four columns! Is your tsv file tab separated?"
					def read2 = row[4] ? row[4].replaceAll(/\]/, "").replaceAll(/ /, "") : null
					def read3 = row[5] ? row[5].replaceAll(/\]/, "").replaceAll(/ /, "") : null
					if (!hasExtension(read1, ".fastq.gz") ) exit 1, "Invalid read names! Reads need to end with .fastq.gz"
					def meta = [:]
					meta.id = row[0]
					def singleEnd = row[1]
					meta.single_end = singleEnd.toBoolean()
					meta.coassemblygroup = row[2]

					if(row[5]){
						return [meta,  [read1, read2, read3] ]
					}else if(row[4]){
						return [meta,  [read1, read2] ]
					}else{
						return [meta,  read1 ]
					}
				}.set { output }

			download_sra(output)
			if(params.getmetadata){
				query_metadata(output)
			}
			reads = download_sra.out.reads
		}else{
			reads = Channel.empty()
		}
	emit:
		reads
}

workflow input_sra_list {
	main: 
		Channel.from(file(params.sralist))
			.splitCsv ( header:false, sep:',' )
			.map {row -> 
					def read1 = row[3] ? row[3].replaceAll(/\[/, "").replaceAll(/\]/, "").replaceAll(/ /, "") : null
					if (!read1) exit 1, "Invalid input samplesheet in at least four columns! Is your tsv file tab separated?"
					def read2 = row[4] ? row[4].replaceAll(/\]/, "").replaceAll(/ /, "") : null
					def read3 = row[5] ? row[5].replaceAll(/\]/, "").replaceAll(/ /, "") : null
					if (!hasExtension(read1, ".fastq.gz") ) exit 1, "Invalid read names! Reads need to end with .fastq.gz"
					def meta = [:]
					meta.id = row[0]
					def singleEnd = row[1]
					meta.single_end = singleEnd.toBoolean()
					meta.coassemblygroup = row[2]

					if(row[5]){
						return [meta,  [read1, read2, read3] ]
					}else if(row[4]){
						return [meta,  [read1, read2] ]
					}else{
						return [meta,  read1 ]
					}
				}.set { output }

		download_sra(output)
		if(params.getmetadata){
			query_metadata(output)
		}
		reads = download_sra.out.reads
	emit:
		reads
}

def hasExtension(it, extension) {
	it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def hasinternetprefix(it) {
	inp = it.toString().toLowerCase()
	inp.startsWith("ftp") || inp.startsWith("http") || inp.startsWith("https") || inp.startsWith("sftp") 
}

class FileCheck {
    def static checkoutfile(def filePath) {
        def file = new File(filePath)
        //Check that file exists and is not empty.
        if (file.exists() && file.isFile() && file.size() > 0) {
            return true  
        } else {
            return false 
        }
    }
}