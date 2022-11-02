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
Reads:          : ${params.reads}
Host genome:    : ${params.genome}"""
if(params.kraken){
log.info "Kraken DB:      : ${params.kraken2_db}"}
if(params.humann){
log.info "HUMAnN DB:      : ${params.humann_db}"}
if(params.humann){
log.info "MPA for HUMAnN  : ${params.metaphlan_db}"}
if(params.metaphlan){
log.info "Metaphlan DB:   : ${params.metaphlan_db}"}
if(params.assembly){
log.info "GTDBTK ref.     : ${params.gtdbtk_reference}"}
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
  nextflow run ikmb/metagenomic-workflows --reads '/path/to/*_R{1,2}_001.fastq.gz' 
  Mandatory arguments:
  --reads 		The path to the fastq.gz files containing PE metagenomic reads (2 per sample) or SE metagenomic reads (1 per sample, use --single_end) or a csv file with a list of reads


  Analysis Modules:
  --metaphlan   Run Metaphlan3 for profiling the composition of microbial communities
                Metaphlan arguments:
                --metaphlan_db  Set directory of metaphlan database (not needed for medcluster)

  --humann      Run HUMAnN3 for profiling the presence/absence and abundance of microbial pathways
                HUMAnN arguments:
                --humann_db     Set directory of humann database (not needed for medcluster)

  --virus       Run Kraken2 with a virus database
                Kraken2 arguments:
                --kraken2_db    Set directory of virus/kraken2 database (not needed for medcluster)

  Experimental:
  --assembly    Run a basic Genome Assembly with Megahit, Metabat and GTDBTK
  --magscot     Run an extended Genome Assembly with 4 binners and MAGScoT for Bin refinement

  Optonal arguments:
  --genome		Remove host contaminations. Use a pre-configured genome sequence by its common name (on medcluster: human, mouse or chimp)
  --cleanreads  Publish QCed fastq.gz files. Disabled by default
  -profile      The nextflow execution profile to use (local or medcluster [default])
  --single_end  Run the pipeline for single-end metagenomic reads

  """.stripIndent()
}
//--email 		An eMail adress to which reports are sent
//TODO: implement: --figures 	Create overview graphics from the result (default: false). Only recommended for smaller sample sizes. 

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

include { MW } from './workflows/metagenomic_workflows' 

//params(params)

workflow {

  MW()

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}