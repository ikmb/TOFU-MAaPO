# Usage:

This pipeline requires Nextflow 21.04.0 or higher. Other dependencies are containerized with Singularity and Docker.<br />

As default, this pipeline uses the profile for Kiel medcluster. Should you choose to run it locally on your own computer, please set **-profile local**. 
Important: Change parameters in conf/local.config to your local hardware specifications prior running the pipeline.


Reference databases for Metaphlan4, Kraken2 and HUMAnN3 (Set also a Metaphlan4 DB for HUMAnN3.6) are needed. On Kiel Medcluster, these are already set in the respective config file.<br />
Metaphlan DB: `--metaphlan_db`<br />
HUMAnN DB:    `--humann_db`<br />
Kraken DB:    `--kraken2_db`<br />
Salmon DB:    `--salmon_db`<br />
Sylph DB:     `--sylph_db`<br />

Pipeline is module based and will run in the most basic run only the QC module.

Run the Pipeline with<br />
```bash
nextflow run ikmb/TOFU-MAaPO --reads '/path/to/fastqfiles/*_R{1,2}_001.fastq.gz'
```
## Input:
Either use:<br />
`--reads` With a glob to your fastq.gz files or to a csv-file containing the columns id, read1,read2 that lists all samples that you want to process. For single-end mode, use only columns "id" and "read1".<br />
or:<br />
`--sra` NCBI SRA Accession ID. Pipeline will download automatically all fastq files for your query. It is mandatory to provide your personal API key for your NCBI account with `--apikey`. Also lists are possible: "--sra ['ERR908507', 'ERR908506', 'ERR908505']". WARNING: The used Nextflow API call to NCBI is not free of bugs. Expect more samples to be processed than are in the input list. Also some samples might be missing. <br />

## Available modules:
For analysis following modules are available:<br />
`--metaphlan` Run Metaphlan4, a tool for profiling the composition of microbial communities<br />
`--humann` Run HUMAnN3, a tool for profiling the abundance of microbial metabolic pathways and other molecular functions<br />
`--kraken` Run Kraken2, a tool for taxonomic classification tool, with a on Medcluster preconfigured RefSeq virus database.<br />
`--bracken` Run Bracken (Bayesian Reestimation of Abundance with KrakEN) after Kraken2. Kraken2 DB must be [bracken-ready](https://github.com/jenniferlu717/Bracken#step-0-build-a-kraken-10-or-kraken-20-database)<br />
`--salmon` Run salmon.<br />
`--sylph` Run sylph.<br />
`--assembly` Run an extended genome assembly workflow with [MAGScoT](https://github.com/ikmb/MAGScoT) Bin Refinement.<br />


## Initialization options:
`--updatemetaphlan` Download the Metaphlan4 database to the directory set in parameter metaphlan_db.<br />
`--updatehumann` Download the HUMAnN3 database to the directory set in parameter humann_db. HUMAnN3 requires the Metaphlan4 database, too.<br />
`--updategtdbtk` Download the GTDB-Tk reference data to the directory set in parameter gtdbtk_reference.<br />


## QC options:
`--genome` Set host genome. On the IKMB Medcluster valid options are human, mouse or chimp. In other cases this needs to be pre-configured. [How to add a host genome to the pipeline?](hostgenome.md) <br />
`--cleanreads`  Publish QC'ed fastq.gz files. Disabled by default.<br /> 
`--no_qc` Skips QC-Module. Only use if your input reads are the output of `--cleanreads`<br /> 
`--fastp` QC is performed with fastp <br /> 

## Metaphlan options:
`--metaphlan_db` Directory of Metaphlan database. REQUIRED! <br /> 
`--publish_metaphlanbam` Publish the bam file output of Metaphlan. <br /> 

## HUMAnN options:
`--metaphlan_db` Directory of Metaphlan database. REQUIRED! <br /> 
`--humann_db` Directory of HUMAnN database. REQUIRED! <br /> 

## Assembly options:
`--assemblymode` Set the mode, if co-assembly (**group** or **all**) or single (**single**, default mode) sample assembly should be performed. The option **group** is only available, if the input is a csv-file with a column "group". In case of co-assembly, only up to 100 samples per group (in "group" mode) or run (in "all" mode) are recommended due to hardware restrictions.<br />
`--binner` Comma separated list of binning tools to use. Options are: **concoct**,**maxbin**,**semibin**,**metabat** and **vamb**. For best performance choose multiple. Default uses all of them. <br />
`--contigsminlength` Set a minimum length of contig. Smaller contigs will be discarded. Default: 2000. <br />
`--semibin_environment` Set the trained environment for SemiBin. Default is **human_gut**. See the [SemiBin Documentation](https://github.com/BigDataBiology/SemiBin/#easy-singleco-assembly-binning-mode) for other options. Choose **global** if no other environment is appropiate.  <br />
`--skip_gtdbtk` Skip GTDB-TK. Both Genome Assembly Modules will run GTDB-TK for taxonomical profiling as a default. <br />
`--skip_checkm` Skip Checkm bin quality check. <br />
`--gtdbtk_reference` GTDB-TK Reference. Reference database for GTDB-TK needs to be set (already set on Kiel Medcluster):<br />
`--publish_megahit` Publish results of megahit with .<br />
`--publish_rawbins` Publish the results of all used binning tools in the genome assembly workflow.<br />
`--vamb_groupsize` Only used when binning with vamb is performed and assemblymode is "single". Set a subgrouping size for vamb, default is 100. This is a temporary fix to enable the pipeline to handle very large cohorts on medium sized hardware. For best results adjust the groupsize to the total sample size of your cohort.<br />
### MAGScoT options:
`--magscot_min_sharing` Scoring parameter a [default=1] <br />
`--magscot_score_a` Scoring parameter a [default=1] <br />
`--magscot_score_b` Scoring parameter b [default=0.5] <br />
`--magscot_score_c` Scoring parameter c [default=0.5] <br />
`--magscot_threshold` Scoring minimum completeness threshold [default=0.5] <br />
`--magscot_min_markers` Minimum number of unique markers in bins to be considered as seed for bin merging [default=25] <br />
`--magscot_iterations` Number of merging iterations to perform. [default=2] <br />

## Other options:
`--single_end` Set the pipeline for single end reads.<br />
`--outdir` Set a custom output directory, default is "results".<br />
`-resume` Resumes pipeline and will continue the run with already completed, cached processes.<br />
`-profile` Change the configuration of the pipeline. Valid options are medcluster (default), local or custom. You can add a new profile for your compute system by editing the file custom.config in the folder conf or create a new one and add it in the file nextflow.config under 'profiles'.<br />
`-work-dir` Set a custom work directory, default is "work".<br />
`-r` Use a specific branch or release version of the pipeline.<br />
`--publish_rawreads` Publish unprocessed/raw files downloaded from SRA in the output directory.<br />

### Kraken2 options:
`--kraken2_db` Directory of used Kraken2 database. Should be Bracken ready for use with Bracken. REQUIRED! <br />

### Salmon options:
`--salmon_db` Directory of used salmon database. REQUIRED! <br />
`--salmon_reference` Path to tab-separated taxonomy file corresponding to the used salmon database. Not required if used with default database. Two column file with header line containing in the first column the bin names used in the salmon database and in the second column the taxonomic assignment by GTDB-Tk in the format "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli". <br />
`--salmon_processing` NOT RECOMMENDED! Shortcut for high-throughput data processing with salmon, skips qc, no other modules available in this mode.  <br />

### Sylph options:
`--sylph_db` Set the path to a sylph databse.<br />
`--sylph_merge` All sylph profiling will be done in one process. Produces a single output for all samples combined. <br />

### Bracken options and their default:
`--bracken_length` = 100<br />
`--bracken_level` = "S"<br />
`--bracken_threshold` = 0<br />
