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
nextflow run ikmb/TOFU-MAaPO --reads '/path/to/fastqfiles/*_R{1,2}_001.fastq.gz'
```
## Available modules:
For analysis following modules are available:<br />
`--metaphlan` Run Metaphlan4, a tool for profiling the composition of microbial communities<br />
`--humann` Run HUMAnN3, a tool for profiling the abundance of microbial metabolic pathways and other molecular functions<br />
`--kraken` Run Kraken2, a tool for taxonomic classification tool, with a on Medcluster preconfigured RefSeq virus database.<br />
`--bracken` Run Bracken (Bayesian Reestimation of Abundance with KrakEN) after Kraken2. Kraken2 DB must be [bracken-ready](https://github.com/jenniferlu717/Bracken#step-0-build-a-kraken-10-or-kraken-20-database)<br />
`--assembly` Run a basic genome assembly workflow.<br />
`--magscot` Run an extended genome assembly workflow with [MAGScoT](https://github.com/ikmb/MAGScoT) Bin Refinement.<br />


## Initialization options:
`--updatemetaphlan` Download the Metaphlan4 database to the directory set in parameter metaphlan_db.<br />
`--updatehumann` Download the HUMAnN3 database to the directory set in parameter humann_db. HUMAnN3 requires the Metaphlan4 database, too.<br />
`--updategtdbtk` Download the GTDB-Tk reference data to the directory set in parameter gtdbtk_reference.<br />


## QC options:
`--genome` Set host genome. On the IKMB Medcluster valid options are human, mouse or chimp. In other cases this needs to be pre-configured.<br />
`--cleanreads`  Publish QC'ed fastq.gz files. Disabled by default.<br /> 
`--no_qc` Skips QC-Module. Only use if your input reads are the output of `--cleanreads`<br /> 

## Metaphlan options:
`--metaphlan_db` Directory of Metaphlan database. REQUIRED!
## HUMAnN options:
`--metaphlan_db` Directory of Metaphlan database. REQUIRED! <br /> 
`--humann_db` Directory of HUMAnN database. REQUIRED! <br /> 
## Assembly options:
`--contigsminlength` Set a minimum length of contig. Smaller contigs will be discarded. Default: 2000. <br />
`--skip_gtdbtk` Skip GTDB-TK. Both Genome Assembly Modules will run GTDB-TK for taxonomical profiling as a default. <br />
`--gtdbtk_reference` GTDB-TK Reference. Reference database for GTDB-TK needs to be set (already set on Kiel Medcluster):<br />
`--publish_megahit` Publish results of megahit with .<br />
`--publish_rawbins` Publish the individual results of all binning tools in the extended genome assembly workflow with.<br />
`--skip_vamb` Deactivate vamb in the magscot workflow. <br />
`--vamb_groupsize` Set a subgrouping size for vamb, default is 100. This is a temporary fix to enable the pipeline to handle very large cohorts on medium sized nodes. For best results adjust the groupsize to the total sample size of your cohort.<br />

## Other options:
`--single_end` Set the pipeline for single end reads.<br />
`--outdir` Set a custom output directory, default is "results".<br />
`-resume` Resumes pipeline and will continue the run with already completed, cached processes.<br />
`-profile` Change the configuration of the pipeline. Valid options are medcluster (default), local or custom. You can add a new profile for your compute system by editing the file custom.config in the folder conf or create a new one and add it in the file nextflow.config under 'profiles'.<br />
`-work-dir` Set a custom work directory, default is "work"<br />
`-r` Use a specific branch or release version of the pipeline.<br />

### Kraken2 options:
`--kraken2_db` Directory of used Kraken2 database. Should be Bracken ready for use with Bracken. REQUIRED! <br />

### Bracken options and their default:
`--bracken_length` = 100<br />
`--bracken_level` = "S"<br />
`--bracken_threshold` = 0<br />
