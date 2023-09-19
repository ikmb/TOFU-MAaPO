/*
 * TOFU-MAaPO
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
TOFU-MAaPO v${workflow.manifest.version}
==========================================
Nextflow Version: $workflow.nextflow.version
Container Engine: ${workflow.containerEngine}
=======INPUTS============================="""
if(params.reads){
log.info "Reads:          : ${params.reads}"}
if(params.sra){
log.info "Accession ID:   : ${params.sra}"}
log.info "Host genome:    : ${params.genome}"
if(params.kraken || params.bracken){
log.info "Kraken DB:      : ${params.kraken2_db}"}
if(params.salmon){
log.info "Salmon DB:      : ${params.salmon_db}"}
if(params.humann){
log.info "HUMAnN DB:      : ${params.humann_db}"}
if(params.metaphlan){
log.info "Metaphlan DB:   : ${params.metaphlan_db}"}
if(params.humann){
log.info "MPdb for HUMAnN : ${params.metaphlan_db}"}
if(params.assembly || params.magscot){
log.info "Assembly mode:  : ${params.assemblymode}"
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
	IKMB | TOFU-MAaPO | v${workflow.manifest.version}
	=================================================================
	Usage:
	The typical command for running the pipeline is as follows:
	nextflow run ikmb/TOFU-MAaPO --reads '/path/to/*_R{1,2}_001.fastq.gz' 
	Mandatory arguments:
	--reads 		The path to the fastq.gz files containing PE metagenomic reads (2 per sample) or SE metagenomic reads (1 per sample, use --single_end) or a csv file with a list of reads
	or
	--sra       Use SRA Accession ID to automatically download the queried fastq files.
				Mandatory:
				--apikey    Your personal NCBI API key.
				Optonial:
				--publish_rawreads	If downloaded SRA files should be copied to the result directry.


	Analysis Modules:
	--metaphlan   Run Metaphlan4 for profiling the composition of microbial communities
				Metaphlan arguments:
				--metaphlan_db    Set directory of metaphlan database (not needed for medcluster)

	--humann      Run HUMAnN3 for profiling the presence/absence and abundance of microbial pathways
				HUMAnN arguments:
				--humann_db       Set directory of humann database (not needed for medcluster)

	--kraken      Run Kraken2
				Kraken2 arguments:
				--kraken2_db      Set directory of virus/kraken2 database (not needed for medcluster)

	--salmon      Run Salmon
				Salmon arguments:
				--salmon_db       Set directory of salmon database

	--assembly    Run a basic Genome Assembly with Megahit, Metabat and GTDB-Tk
	--magscot     Run an extended Genome Assembly with 4 binning tools and MAGScoT for bin refinement.
				Assembly arguments:
				--gtdbtk_reference Set directory of the GTDB-Tk reference data

	Optional arguments:
	General:
		-profile                The nextflow execution profile to use (custom, local or medcluster [default])
		--single_end            Run the pipeline for single-end metagenomic reads
		-work-dir               Set a custom work directory, default is "work".
		--outdir                Set a custom work directory for all outputs, default is "results".
		-resume                 Resumes pipeline and will continue the run with already completed, cached processes.

	Initialization:
		--updatemetaphlan       Download the Metaphlan4 database to the directory set in parameter metaphlan_db
		--updatehumann          Download the HUMAnN3 database to the directory set in parameter humann_db. HUMAnN3 requires the Metaphlan4 database, too.
		--updategtdbtk          Download the GTDB-Tk reference data to the directory set in parameter gtdbtk_reference.

	QC:
		--genome		            Remove host contaminations. Use a pre-configured genome sequence by its common name (on medcluster: human, mouse or chimp)
		--cleanreads            Publish QCed fastq.gz files. Disabled by default.
		--no_qc                 Skips QC-Module. Only use if your input reads are the output of --cleanreads. Not recommended.

	Kraken2/Bracken:
		--kraken2_db            Directory of used Kraken2 database. Should be Bracken ready for use with Bracken. REQUIRED!
		--bracken_length        Read length. Default: 100
		--bracken_level         Taxonomic level. Options: D,P,C,O,F,G,S,S1. Default: "S"
		--bracken_threshold     Number of reads required prior to abundance estimation. Default: 0
	
	Metaphlan:
		--metaphlan_db          Directory of Metaphlan database. REQUIRED!

	HUMAnN:
		--humann_db             Directory of HUMAnN database. REQUIRED!
		--metaphlan_db          Directory of Metaphlan database. REQUIRED!

	Assembly:
		--assemblymode			Select an assembly mode, can be "single" or for co-assembly "all". Also possible: "group" if input is a csv-file that contains a column "group"
		--contigsminlength	    Minimum length of contigs. Default: 2000
		--gtdbtk_reference      Directory of database. REQUIRED for GTDB-Tk!
		--skip_gtdbtk           Skip GTDB-Tk.
		--publish_megahit       Publish all megahit contigs. Default: false
		--publish_rawbins       Publish all raw bins from all binning tools. Default: false
		--contig_sep            Contig name separator. Default: "megahitcontig"
		--skip_vamb             Skip Vamb.
		--vamb_groupsize        Group size of samples to use in one Vamb run if assembly mode is "single". Recommended: Number of all samples in the run. Default: 100.
	""".stripIndent()
}
//--email 		An eMail adress to which reports are sent
//TODO: implement: --figures 	Create overview graphics from the result (default: false). Only recommended for smaller sample sizes. 

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

include { tofumaapo } from './workflows/tofumaapo' 
include { tofumaapo_salmon } from './workflows/tofumaapo_salmon' 
//params(params)

workflow {
	if(!params.salmon_processing){
		tofumaapo()
	}else{
		tofumaapo_salmon()
	}

}

workflow.onComplete {
	log.info "========================================="
	log.info "Duration:		$workflow.duration"
	log.info "========================================="
}