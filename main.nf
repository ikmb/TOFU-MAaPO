/*
 * TOFU-MAaPO
 * By Eike Matthias Wacker (e.wacker@ikmb.uni-kiel.de)
 * Based on https://github.com/marchoeppner/metagenomic-profiling Release 1.2 by Marc P. Hoeppner (m.hoeppner@ikmb.uni-kiel.de)
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2



include { tofumaapo } from './workflows/tofumaapo' 
include { tofumaapo_sylph } from './workflows/tofumaapo_sylph' 
//params(params)

workflow {
	main:
	WorkflowHeaderHelp.initialise(workflow, params, log)
	WorkflowCheck.startuptests(workflow, params, log)

	if(params.sylph_processing){
		tofumaapo_sylph()
	}else{
		tofumaapo()
	}

	workflow.onComplete = {
	log.info "========================================="
	log.info "Pipeline startet at:	$workflow.start"
	log.info "Pipeline completed at:	$workflow.complete"
	log.info "Duration:		$workflow.duration"
	log.info "Execution status:	${ workflow.success ? 'OK' : 'failed' }"
	log.info "========================================="
	}
}