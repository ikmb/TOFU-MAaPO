# Installation
Before you can run TOFU-MAaPO, you need to install [Singularity](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html) and [Nextflow](https://www.nextflow.io/docs/latest/install.html). We show how to install the dependencies in the [the Quick start section](../README.md#installing-dependencies) or you install them manually by following the links.<br />

You will need to prepare databases for respective modules (e.g., Metaphlan4, HUMAnN3) and a config file for your compute system.<br />
# Configuration

Before running **TOFU-MAaPO** for the first time, you must set up a configuration file tailored to your computing system. There are two approaches for configuration:

### Option 1: Modify the repository's config file (not recommended)
- Edit the [`custom.config`](../conf/custom.config) file located in the repository's `conf` directory.
- Run the pipeline using the custom profile:
  ```bash
  -profile custom
  ```
>**Caution**: This approach is less flexible and may cause issues when updating the repository, as your changes might be overwritten.

### Option 2: Create a separate configuration file (recommended)
- Create and maintain a separate configuration file (e.g., `tofu.config`).
- Include your custom configuration in the pipeline execution by adding:
  ```bash
  -profile custom -c /path/to/tofu.config
  ```

In this example `tofu.config` should contain following (from you customized) lines: <br />
```
params {
  //Software DB locations, UNCOMMENT AND CHANGE THEM:

  //metaphlan_db = "${baseDir}/databases/Metaphlan/4.0"
  //kraken2_db = "${baseDir}/databases/Kraken2/k2_viral_20210517"
  //humann_db = "${baseDir}/databases/Humann3/3.6" 
  //gtdbtk_reference = "${baseDir}/databases/GTDB-TK/release207_v2"

  //For host read removal list your host genomes as bowtie2 index in this named list with full path to the basename of the index:
	'genomes' {
		'human' {
      bowtie_index = false
		}
	}

  // MAXIMUM PER PROCESS CONFIGS, CHANGE THEM TO YOUR HARDWARE SPECS
  max_memory = 100.GB
  max_cpus = 32
  max_time = 48.h
  
  //Scratch. Does your system support scratch? If not, leave it false
  scratch = false
}

//Enable Singularity as container software
singularity {
	enabled = true
  // Singularity configs, CHANGE THEM TO YOUR USED FILESYSTEM, if not properly set the container won't see your files
	runOptions = "-B /home -B /tmp"
	// where should the containers be downloaded to
	cacheDir = "${launchDir}/singularity_cache"
}

process {
  //Default executor for each process, other options can be e.g. SLURM.
  //See https://www.nextflow.io/docs/latest/executor.html for more options and details.
  executor='local'
}

//Default for total execution, remove this whole part if you are using a different option above than executor='local':
executor {
  cpus = 32
  memory = 200.GB
}
```

### Executors

If you intend to run the pipeline on an HPC or cloud service, refer to the [Nextflow Executors documentation](https://www.nextflow.io/docs/latest/executor.html) for guidance on adapting your custom configuration file to your environment.

---

### QC: Host Decontamination

To perform host decontamination, you need the appropriate host genome Bowtie2 indexes. Download these indexes, for example, from [Bowtie2 AWS indexes](https://benlangmead.github.io/aws-indexes/bowtie), and specify the path to the basename of the index files in your custom configuration file before running the pipeline. Add the following snippet to your configuration file:

```groovy
'genomes' {
    'human' {
        bowtie_index = "/path/to/your/references/iGenomes/references/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome"
    }
}
```
### MetaPhlAn & HUMAnN
The pipeline requires databases to run HUMAnN and MetaPhlAn.

Use the pipeline's built-in download functionality for HUMAnN and MetaPhlAn databases with the parameters:

`--updatemetaphlan` for MetaPhlAn4.  
`--updatehumann` for HUMAnN3.  

Supply the paths to store the databases to using:  
`--metaphlan_db` for MetaPhlAn4.  
`--humann_db` for HUMAnN3.

Run the pipeline with these parameters **during your initial run**.   
Ensure your local computer has internet access, as the pipeline will automatically download and extract the databases. 
> Hint: Once the first run is successful, update your custom configuration file to permanently set the paths for `--metaphlan_db` and `--humann_db`.

>**Note**: Running HUMAnN3 also requires the MetaPhlAn4 database (version vJan21).

Update your configuration file with the following snippet:
```groovy
params {
    metaphlan_db = "/path/to/your/databases/Metaphlan/4.0"
    humann_db = "/path/to/your/databases/Humann3/3.6"
}
```

### Kraken2
Download and extract a Kraken2 database, such as those available from [Kraken2 AWS indexes](https://benlangmead.github.io/aws-indexes/k2). You can specify the database path either in your custom configuration file or dynamically during each run with the parameter `--kraken2_db`.

>**Note**: If you plan to use **Bracken**, ensure the Kraken2 database is Bracken-ready.

Add the following snippet to your custom configuration file:
```groovy
params {
    kraken2_db = "/path/to/your/databases/Kraken2/krakendb"
}
```

### Genome Assembly
For genome assembly, the pipeline requires GTDB-Tk reference data (currently Release R207_v2). 
Use the pipeline's built-in download functionality to download the GTDBTk database with the parameter 
`--updategtdbtk` and specify the path where to store the database with `--gtdbtk_reference`.

Add the following snippet to your custom configuration file:
```groovy
params {
	gtdbtk_reference = "/path/to/your/databases/GTDB-TK/release207_v2"
}
```

# Kiel CAUcluster
On Kiel Caucluster, please load the following modules with:
```bash
module load gcc12-env
module load singularity nextflow
```
When starting TOFU-MAaPO, make sure the CAUcluster profile is selected, by using: `-profile caucluster`.
No further action or initialization steps are required for IKMB users! <br />
Databases for all tools (with the exception of Salmon) and the host genomes human, mouse and chimp are already set up.
>Note: Due to permission-restrictions all non-IKMB users need to initialize the pipeline and download their databases in their own directories. Create a configuration file as described above for all database and reference paths. You do not need to set the other configurations, but can simply use `-profile caucluster`.  <br />

> Note: It is recommended to start the pipeline in a tmux or screen session directly on one of the head nodes. This avoids SLURM timelimits, which could kill your pipeline run before it finishes, if you start the pipeline via a SLURM script. <br />