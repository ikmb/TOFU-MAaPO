/*
 * METAGENOMIC W O R K F L O W S
 * By Eike Matthias Wacker (e.wacker@ikmb.uni-kiel.de)
 * Based on https://github.com/marchoeppner/metagenomic-profiling Release 1.2 by Marc P. Hoeppner (m.hoeppner@ikmb.uni-kiel.de)
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters in nextflow.config
 */ 
log.info """\
METAGENOMIC W O R K F L O W S  v${workflow.manifest.version}
==========================================
Nextflow Version: $workflow.nextflow.version
Container Engine: ${workflow.containerEngine}
=======INPUTS=============================
Reads:          : ${params.reads}"
Metaphlan DB:   : ${params.metaphlan_db}
HUMAnN DB:      : ${params.humann_db}
Kraken DB:      : ${params.kraken2_db}
Host genome:    : ${params.genome}"""
log.info "=========================================="
log.info "Command Line:     $workflow.commandLine"
log.info "=========================================="

/*
 * Show help message with "--help"
 */
def helpMessage() {
  log.info"""
  =================================================================
   IKMB | METAGENOMIC W O R K F L O W S | v${workflow.manifest.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run eikematthias/metagenomic-workflows --reads '/path/to/*_R{1,2}_001.fastq.gz' 
  Mandatory arguments:
  --reads 		The path to the fastq.gz files containing PE metagenomic reads (1 per sample)

  Analysis Modules:
  --metaphlan   Run Metaphlan3 for profiling the composition of microbial communities
                Metaphlan arguments:
                --metaphlan_db  Set directory of metaphlan database (not needed for medcluster)

  --humann      Run HUMAnN3 for profiling the presence/absence and abundance of microbial pathways
                HUMAnN arguments:
                --humann_db     Set directory of humann database (not needed for medcluster)

  --virus       Run Kraken2 with a virus database
                Kraken2 argmeunts:
                --kraken2_db    Set directory of virus database (not needed for medcluster)

  Experimental:
  --assembly    Run Genome Assembly with Megahit, Metabat and GTDBTK

  Optonal arguments:
  --genome		Remove host contaminations. Use a pre-configured genome sequence by its common name (on medcluster: human, mouse or chimp)
  --email 		An eMail adress to which reports are sent
  -profile      The nextflow execution profile to use (local or medcluster [default])

  """.stripIndent()
}
//TODO: implement: --figures 	Create overview graphics from the result (default: false). Only recommended for smaller sample sizes. 

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

/* 
 * Import QC modules
 * MULTIQC1 and MULTIQC2 differ only in having 1 or 2 input-channels, remove this when optional inputs are possible
 */
include { 
  FASTQC;
  TRIMREADS;
  CLEANPEREADS;
  CLEANSEREADS;
  FILTERREADS;
  MULTIQC1;
  MULTIQC2
  } from './modules/QC.nf' 

/* 
 * QC pipeline logic
 * QC_rmHost removes all Host contaminations, QC_noHost does not filter for Host contaminations
 */
workflow QC_rmHost{
    main:
        if (!params.genome) {
            println "No Host genome was specified!"
        } else if (!params.genomes) {
            exit 1, "Specified a genome name for host mapping, but no genomes are configured for your profile. Aborting!"
        } else if (!params.genomes.containsKey(params.genome)) {
            exit 1, "Specified unknown name for the host genome! Valid options are: ${params.genomes.keySet()}"
        } else {
            log.info "Using ${params.genome} as host species."

            bowtie_base = params.genomes[params.genome].bowtie_index
            //bloom_index = params.genomes[params.genome].bloom_index
        }

     reads_ch = Channel.fromFilePairs(params.reads , flat: true ).ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
     FASTQC(reads_ch)
     TRIMREADS(reads_ch)
     CLEANPEREADS(TRIMREADS.out.filterPEReads)
     CLEANSEREADS(TRIMREADS.out.filterSEReads)
     FILTERREADS(
         CLEANPEREADS.out.join(CLEANSEREADS.out),
         Channel.fromPath("${bowtie_base}*").collect(),         
         Channel.fromPath(params.genomes[params.genome].bowtie_index).map{index -> index.Name} )    
    emit:
     my_data = FILTERREADS.out.cleaned_reads
     fastqcoutput = FASTQC.out.collect()
}

workflow QC_noHost{
    main:
        
    println "No Host genome was specified!"

     reads_ch = Channel.fromFilePairs(params.reads , flat: true ).ifEmpty {exit 1, "Could not find the specified input reads $params.reads"}
     FASTQC(reads_ch)
     TRIMREADS(reads_ch)
     CLEANPEREADS(TRIMREADS.out.filterPEReads)
     CLEANSEREADS(TRIMREADS.out.filterSEReads)
     
	emit:
     my_data = CLEANPEREADS.out.join(CLEANSEREADS.out)
     fastqcoutput = FASTQC.out.collect()    
}
/* 
 * Import kraken modules 
 */
  include {
  KRAKEN2;
  KRAKEN2YAML
  } from './modules/kraken.nf' 
  
/*
 * Kraken2 pipeline logic
 */
workflow kraken{
    take: data
    main:
        KRAKEN2(data)
        KRAKEN2YAML(KRAKEN2.out.krakenreport.collect())
    emit:
        kraken_data = KRAKEN2YAML.out
}

/* 
 * Import metaphlan modules 
 */
  include {
  PREPARE_METAPHLAN;
  METAPHLAN;
  ABUNDANCEMERGE
  } from './modules/metaphlan.nf'

/*
 * Metaphlan3 pipeline logic
 */
workflow metaphlan{
    take: data
    main:
        /*
        if (params.metaphlan_db) {
            db_path = file("${params.metaphlan_db_test}")
            if (!db_path.exists()) { exit 1, "Could not find Metaphlan database - please check the path" } 
        }else {exit 1, "No Metaphlan database was specified, aborting..."}
        */
        PREPARE_METAPHLAN()
        METAPHLAN(data)
        ABUNDANCEMERGE(METAPHLAN.out.outputMetaphlan.collect() )
}

/*
 * Import humann modules
 */
  include {
      HUMANN;
      JOINgenefamilies;
      JOINpathabundance;
      JOINpathcoverage
  } from './modules/humann.nf'

/*
 * Humann3 pipeline logic
 */
workflow humann{
    take: data
    main:
        /*
        if (params.metaphlan_db) {
            db_path = file("${params.metaphlan_db_test}")
            if (!db_path.exists()){ exit 1, "Could not find Metaphlan database - please check the path"}
        } else {
	    exit 1, "No Metaphlan database was specified, aborting..."
        }
        //TODO: CHANGE THIS for a path checking solution:
        
        humann_db_path = file("${params.humann_db_test}")
	    if (!humann_db_path.exists() ){ exit 1, "Could not find your HUMAnN DB - please check the path"
        } else if (params.humann) {
        exit 1, "No HUMAnN database was specified, aborting..."
        }
        */
        if(!params.metaphlan){
            PREPARE_METAPHLAN()
        }
        HUMANN(data)
        JOINgenefamilies(HUMANN.out.genefamilies.collect() )
        JOINpathabundance(HUMANN.out.pathabundance.collect())
        JOINpathcoverage(HUMANN.out.pathcoverage.collect())
}

/*
 * Import Genome Assembly modules
 */
 include {
     MEGAHIT;
     MAPPING;
     METABAT;
     contigs_to_bins;
     checkm_all_bins;
     GTDBTK;
     getCountTable
 } from './modules/assembly.nf'

/*
 * Genome Assembly pipeline logic
 */
workflow assembly{
    take: data
    main:
        MEGAHIT(data)
        MAPPING(MEGAHIT.out.contigs)
        METABAT(MAPPING.out.maps)
        contigs_to_bins(METABAT.out)
        checkm_all_bins(METABAT.out)
        GTDBTK(METABAT.out)
        getCountTable(MAPPING.out.counttable)
}

/* 
 * Main pipeline logic
 */
workflow {
    main:
    //QC, either with host contaminations or without:
        if(params.genome){
            QC_rmHost()
            QCout = QC_rmHost.out.my_data
            Fastqcoutput = QC_rmHost.out.fastqcoutput.collect()
        }else{
            QC_nohost()
            QCout = QC_noHost.out.my_data
            Fastqcoutput = QC_noHost.out.fastqcoutput.collect()
            }
    //kraken:
        if(params.virus){
        kraken(QCout)
        }
    //metaphlan:
        if(params.metaphlan){
            metaphlan(QCout)
        }
    //humann:
        if(params.humann){
            humann(QCout)
        }
    //genome assembly:
        if(params.assembly){
            assembly(QCout)
        }
    //multiqc, collecting all fastqc- and kraken-files; change this when optional inputs are doable:
        if(params.virus){
            MULTIQC2(
                Fastqcoutput.collect(),
                kraken.out.collect()
            )    
        }else{
            MULTIQC1(
                Fastqcoutput.collect()
            ) 
        }
}