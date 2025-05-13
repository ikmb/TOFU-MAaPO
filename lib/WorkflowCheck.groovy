
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
                if(!FileCheck.isValid("${params.sylph_db}") ){
                    log.info "The sylph database is either not existing or is corrupted please check that the path to the database '${params.sylph_db}' is valid, exiting now."
                    System.exit(1)
                }
            }else{
                log.info "No sylph database directory was specified, please specify a sylph database directory with '--sylph_db'. Exiting now."
                System.exit(1)
            }
        }

        if(params.updatemetaphlan){
            if(!params.metaphlan_db ){
                    log.info "The location for the metaphlan database is not set! Set a path for the location with '--metaphlan_db', exiting now."
                    System.exit(1)
                }
        }

        if(params.updatehumann){
            if(!params.humann_db ){
                    log.info "The location for the humann database is not set! Set a path for the location with '--humann_db', exiting now."
                    System.exit(1)
                }
        }

        if(params.genome){
            if(!FileCheck.checkoutmultFiles("${params.genomes[params.genome].bowtie_index}") ){
                log.info "The selected host genome '${params.genome}' with the prefix '${params.genomes[params.genome].bowtie_index}' cannot be found. Is the path and basename correct?"
                System.exit(1)
            }
        }
    }

}

import java.net.URL
import java.net.HttpURLConnection
import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.Path
import java.nio.file.Files
import java.nio.file.DirectoryStream
import java.io.IOException
class FileCheck {
    /*
        Static method to validate if the input string is a locally available, valid file.
        
        @param filePath A string that looks like a file path.
        @return boolean True if the input is a valid file (exists, is a file, and is > 0KB).
     */
    def static checkoutfile(def filePath) {
        def file = new File(filePath)
        //Check that file exists and is not empty.
        if (file.exists() && file.isFile() && file.size() > 0) {
            return true  
        } else {
            return false 
        }
    }

    /*
        Static method to validate if the input string is a valid file or a valid URL.
        
        param input A string that looks like a file path or a URL (expected is online available database).
        return boolean True if the input is a valid file (exists, is a file, and is > 0KB) or a valid URL.
     */
    static boolean isValid(String input) {
        if (isFile(input)) {
            return isFileValid(input)
        } else if (isURL(input)) {
            return isURLValid(input)
        }
        return false
    }

    /*
        Check if the input string is a valid file.
        
        param filePath The file path to check.
        return boolean True if the file exists, is a file, and is greater than 0KB.
     */
    private static boolean isFileValid(String filePath) {
        def path = Paths.get(filePath)
        return Files.exists(path) && Files.isRegularFile(path) && Files.size(path) > 0
    }

    /*
        Check if the input string is a URL.
        
        param input The input string to check.
        return boolean True if the input is a valid URL format.
     */
    private static boolean isURL(String input) {
        try {
            new URL(input)
            return true
        } catch (MalformedURLException e) {
            return false
        }
    }

    /*
        Check if the input string is a valid URL.
        
        param urlString The URL string to check.
        return boolean True if the URL is accessible (HTTP response code 200).
     */
    private static boolean isURLValid(String urlString) {
        try {
            URL url = new URL(urlString)
            HttpURLConnection connection = (HttpURLConnection) url.openConnection()
            connection.setRequestMethod("HEAD")
            connection.setConnectTimeout(3500) // 3.5 seconds timeout
            connection.connect()
            int responseCode = connection.getResponseCode()
            return responseCode == HttpURLConnection.HTTP_OK
        } catch (Exception e) {
            return false
        }
    }

    /*
        Check if the input string is a file path.
        
        param input The input string to check.
        return boolean True if the input is likely a file path.
     */
    private static boolean isFile(String input) {
        return Files.exists(Paths.get(input))
    }

    /*
        Check if a given path plus basename corresponds to one or more existing files,
        of which at least one is a file larger than 0KB.
    
        param filePath The input string representing the path plus basename.
        return true if the condition is met, otherwise false.
     */
    static boolean checkoutmultFiles(String filePath) {
        def files = []
        Path pathObj = Paths.get(filePath)
        Path directory = pathObj.getParent() ?: Paths.get(".")  // Use current directory if no parent
        String basename = pathObj.getFileName().toString()
        if (!basename.contains("*")) {
            basename = basename + "*"
        }
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(directory, basename)) {
            for (Path file : stream) {
                files << file
            }
        } catch (IOException e) {
            e.printStackTrace()
            return false
        }

        return files.any { file ->
            Files.exists(file) && Files.isRegularFile(file) && Files.size(file) > 0
        }
    }
}