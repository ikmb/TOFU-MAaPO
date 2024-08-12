# Installation

This pipeline is containerized with Singularity. 

Before you can run TOFU-MAaPO you need to install [Singularity](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html) and [Nextflow](https://www.nextflow.io/docs/latest/install.html). We show how to install the dependencies in the [the Quick start section](../README.md#quick-start) or you install them manually by following the links.<br />

You will need to prepare databases for respective modules (Metaphlan4, HUMAnN3 or Kraken2) and a config file for your compute system.<br />

# Configuration:
Before you run TOFU-MAaPO the first time you need to create a configuration file for your computing system. You can either choose to change the [custom.config file in the repository directory conf](../conf/custom.config). Then you will be able to run the pipeline with your config with the parameter `-profile custom`. A second option is to create a separate file that you can include into your Nextflow call with the addition of `-profile custom -c tofu.config`. <br /><br />

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
Should you want to run the pipeline on a HPC or Cloud service, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html) to adapt your custom config file.

### QC
For host decontamination: Download your needed host genome as Bowtie2 indexes from e.g. [here](https://benlangmead.github.io/aws-indexes/bowtie) and set the the path to the  basename of the index files in your custom config file prior running the pipeline like so:
```
'genomes' {
	'human' { bowtie_index = "/path/to/your/references/iGenomes/references/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome"
	}
}
```
### Metaphlan & Humann
The pipeline needs databases for HUMAnN3, Metaphlan4 and Kraken2. For HUMAnN3 and Metaphlan4 the pipeline can download the necessary files with the parameters `--updatemetaphlan` and `--updatehumann` to the paths you supply with `--metaphlan_db` and `--humann_db`. 
Run the pipeline with the given parameters in your first run. The pipeline will then download and unzip the databases at this path. Make sure your local computer is connected with the internet. After the pipeline sucessfully finished the first run, change your custom config file to set metaphlan_db and humann_db to your now downloaded databases. Note that for running HUMAnN3 besides of the HUMAnN3 database you will also need the Metaphlan4 database (in Version vJan21).

Edit and include this snipped in your custom config:
```
params {
	metaphlan_db = "/path/to/your/databases/Metaphlan/4.0"
	humann_db = "/path/to/your/databases/Humann3/3.6"
}
```

### Kraken2
Download and extract a Kraken 2 database for example from [here](https://benlangmead.github.io/aws-indexes/k2) and add them to your custom config file or set a path each run with the parameter `--kraken2_db`. Make sure, this database is also Bracken ready if you want to use Bracken.

Edit and include this snipped in your custom config:
```
params {
	kraken2_db = "/path/to/your/databases/Kraken2/k2_viral_20210517"
}
```

### Genome assembly
For genome assembly: Please download and extract the GTDB-Tk reference data (currently Release R207_v2) from [here](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data) and set the path either in your custom config or set it with parameter `--gtdbtk_reference`. The needed data can also be automatically installed with `--updategtdbtk` to the path that is supplied by `--gtdbtk_reference`.

Edit and include this snipped in your custom config:
```
params {
	gtdbtk_reference = "/path/to/your/databases/GTDB-TK/release207_v2"
}
```

# Kiel Medcluster
On Kiel Medcluster, please load the following modules with:
```bash
module load singularity nextflow
```
When starting TOFU-MAaPO, make sure the medcluster profile is selected, by using: `-profile medcluster`
No further action is required! <br />
Databases for all tools (with the exception of Salmon) and the host genomes human, mouse and chimp are already set up.