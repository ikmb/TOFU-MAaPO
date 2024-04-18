
class WorkflowCheck {
    
    public static void startuptests(workflow, params, log) {

        if (params.salmon) {
            if(params.salmon_db){
                if(!FileCheck.checkoutfile("${params.salmon_db}/reflengths.bin")){
                    log.info "The salmon database is either not existing or is corrupted! Please check that the path to the database '${params.salmon_db}' is valid. Exiting now."
                    System.exit(1)
                }
            }else{
                log.info "No salmon database directory was specified. Please specify a salmon database directory with '--salmon_db'. Exiting now."
                System.exit(1)
            }
        }

        if(params.sylph){
            if(params.sylph_db ){
                if(!FileCheck.checkoutfile("${params.sylph_db}") ){
                    log.info "The sylph database is either not existing or is corrupted please check that the path to the database '${params.sylph_db}' is valid, exiting now."
                    System.exit(1)
                }
            }else{
                log.info "No sylph database directory was specified, please specify a sylph database directory with '--sylph_db'. Exiting now."
                System.exit(1)
            }
        }

        if (!(params.assemblymode == "single" || params.assemblymode == "all" || params.assemblymode == "group")) {
            log.info  "Parameter 'assemblymode' can only be set to 'single', 'all' or 'group', exiting now."
            System.exit(1)
        }

        if('spades' in params.assembler.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '').replaceAll('2', '').replaceAll('meta', '')}){
            if(params.assemblymode == "all" || params.assemblymode == "group"){
                log.info  "Parameter 'assemblymode' can only be set to 'single', 'all' or 'group', exiting now."
                System.exit(1)
            }

            if(params.coassembly){
                log.info  "Parameter 'coassembly' can not be used together with '--spades', exiting now."
                System.exit(1)
            }
        }

    }

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