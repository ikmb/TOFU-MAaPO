# Usage:

This pipeline requires Nextflow 21.04.0 or higher. Other dependencies are containerized with Singularity and Docker.<br />

As default, this pipeline works with the profile for Kiel medcluster. Should you choose to run it locally on your own computer, please set **-profile local**. 
Important: Change parameters in conf/local.config to your local hardware specifications prior running the pipeline.

On Kiel Medcluster, please load the following modules with:
```bash
module load singularity nextflow
```

Reference databases for Metaphlan4, Kraken2 and HUMAnN3 (Set also aMetaphlan4 DB for HUMAnN3.6) are needed. On Kiel Medcluster, these are already set in the respective config file.<br />
Metaphlan DB: `--metaphlan_db`<br />
HUMAnN DB:    `--humann_db`<br />
Kraken DB:    `--kraken2_db`<br />

Pipeline is module based and will run in the most basic run only the QC module.

Run the Pipeline with<br />
```bash
nextflow run ikmb/metagenomic-workflows --reads '/path/to/fastqfiles/*_R{1,2}_001.fastq.gz'
```
## Available modules:
For analysis following modules are available:<br />
**--metaphlan** Run Metaphlan4, a tool for profiling the composition of microbial communities<br />
**--humann** Run HUMAnN3, a tool for profiling the abundance of microbial metabolic pathways and other molecular functions<br />
**--virus** Run Kraken2, a tool for taxonomic classification tool, with a on Medcluster preconfigured RefSeq virus database.<br />
**--bracken** Run Bracken (Bayesian Reestimation of Abundance with KrakEN) after Kraken2. Kraken2 DB must be [bracken-ready](https://github.com/jenniferlu717/Bracken#step-0-build-a-kraken-10-or-kraken-20-database)<br />

Experimental:<br />
**--assembly** Run a basic genome assembly workflow.<br />
**--magscot** Run an extended genome assembly workflow with [MAGScoT](https://github.com/ikmb/MAGScoT) Bin Refinement.<br />
Both Genome Assembly Modules will run GTDB-TK for taxonomical profiling.<br />
Skip GTDB-TK with `--skip_gtdbtk`<br />
Reference database for GTDB-TK needs to be set (already set on Kiel Medcluster):<br />
GTDB-TK Reference: `--gtdbtk_reference`<br />


## QC options:
**--genome** Set host genome. On the IKMB Medcluster valid options are human, mouse or chimp. In other cases this needs to be pre-configured.<br />
**--cleanreads**  Publish QC'ed fastq.gz files. Disabled by default.<br /> 
**--no_qc** Skips QC-Module. Only use if your input reads are the output of `--cleanreads`<br /> 

## Other options:
**--single_end** Set the pipeline for single end reads.<br />
**--outdir** Set a custom output directory, default is "results".<br />
**-resume** Resumes pipeline and will continue the run with already completed, cached processes.<br />
**-profile** Change the configuration of the pipeline. Valid options are medcluster (default), local or dx. You can add a new profile for your compute system by adding a new config file in the folder conf and adding it in the file nextflow.config under 'profiles'.<br />
**-work-dir** Set a custom work directory, default is "work"<br />
**--updatemetaphlan** Check whether metaphlan-db is still up-to-date before running. Update must be made manually<br />

### Bracken options and their default:
**--bracken_length** = 100<br />
**--bracken_level** = "S"<br />
**--bracken_threshold** = 0<br />
