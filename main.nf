/*
 * TOFU-MAaPO
 * By Eike Matthias Wacker (e.wacker@ikmb.uni-kiel.de)
 * Based on https://github.com/marchoeppner/metagenomic-profiling Release 1.2 by Marc P. Hoeppner (m.hoeppner@ikmb.uni-kiel.de)
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

WorkflowHeaderHelp.initialise(workflow, params, log)
WorkflowCheck.startuptests(workflow, params, log)


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