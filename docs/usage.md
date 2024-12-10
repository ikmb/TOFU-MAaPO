# TOFU-MAaPO Usage Guide

>**Note:** You must create a [**custom configuration file**](installation.md#configuration) and add it with `-profile custom -c tofu.config` to your TOFU-MAaPO call!

## Reference Databases
Reference databases are mandatory for specific modules:  
MetaPhlAn DB: `--metaphlan_db`<br />
HUMAnN DB:    `--humann_db`<br />
Kraken2 DB:    `--kraken2_db`<br />
Sylph DB:     `--sylph_db`<br />
Salmon DB:    `--salmon_db`<br />
GTDB-Tk DB:    `--gtdbtk_reference`<br />

## Initialization options
> **Parameters should be set in the first run of the respective module**  

For **MetaPhlAn**, **HUMAnn** and **GTDB-Tk** the pipeline can download the required database via following flags:  

- `--updatemetaphlan` Download the Metaphlan4 database to the directory set in parameter `--metaphlan_db`.<br />
- `--updatehumann` Download the HUMAnN3 database to the directory set in parameter `--humann_db`. HUMAnN3 requires the Metaphlan4 database, too.<br />
- `--updategtdbtk` Download the GTDB-Tk reference data to the directory set in parameter `--gtdbtk_reference`.<br />

## Basic Execution
By default, only the **Quality Control** runs unless additional modules are specified.  

Example command:
```bash
nextflow run ikmb/TOFU-MAaPO --reads '/path/to/fastqfiles/*_R{1,2}_001.fastq.gz' -profile custom -c tofu.config
```
## Input Options
Choose one of the following options:
### 1. FASTQ Files  
Use the `--reads` parameter with a **glob pattern** for your .fastq.gz files as seen above **or** 
provide a **CSV file** with the following columns:

- `id`: Sample identifier  
- `read1`: Path to the forward reads
- `read2`: Path to the reverse reads (for paired-end data)  
  
For single-end reads, include only `id` and `read1`.

### 2. SRA Accessions
Provide SRA Accession IDs via the `--sra` option.  <br />
**Mandatory: Provide your personal NCBI API key with `--apikey`.**   <br />
The pipeline will automatically download the corresponding FASTQ files. Example:
```bash
--sra 'SRX1234567' --apikey **YOUR_NCBI_API_KEY**
```
For mulitple IDs, use:
```bash
--sra ['ERR908507', 'ERR908506', 'ERR908505'] --apikey **YOUR_NCBI_API_KEY**
```
> **Note**: The Nextflow API call to NCBI may result in extra or missing samples. Ensure to verify downloaded data.

# Available modules

For analysis following modules are available:<br />
## Genome assembly
`--assembly` Run an extended genome assembly workflow with [MAGScoT](https://github.com/ikmb/MAGScoT) bin refinement.<br />

## Assemby-free metabolic gene abundance estimation
`--humann` Run HUMAnN3, a tool for profiling the abundance of microbial metabolic pathways and other molecular functions<br />
## Taxonomical abundance tools
- `--metaphlan` Run **MetaPhlAn4**, a tool for profiling the composition of microbial communities<br />
- `--kraken` Run **Kraken2**, a tool for taxonomic classification tool.<br />
- `--bracken` Run **Bracken** (Bayesian Reestimation of Abundance with KrakEN) after Kraken2. Kraken2 DB must be [bracken-ready](https://github.com/jenniferlu717/Bracken#step-0-build-a-kraken-10-or-kraken-20-database)<br />
- `--sylph` Run **Sylph**.<br />
- `--salmon` Run **Salmon**. *Usage not recommended* <br />



# General options
`--outdir` Set a custom output directory, default is "results".<br />
`-resume` Resumes pipeline and will continue the run with already completed, cached processes.<br />
`-profile` Change the configuration of the pipeline. Valid options are medcluster (default), local or custom. You can add a new profile for your compute system by editing the file custom.config in the folder conf or create a new one and add it in the file nextflow.config under 'profiles'.<br />
`-work-dir` Set a custom work directory, default is "work".<br />
`-r` Use a specific branch or release version of the pipeline.<br />
`--publish_rawreads` Publish unprocessed/raw files downloaded from SRA in the output directory.<br />

# Module specific options
## QC options
- `--cleanreads`  Save QC'ed FASTQ files (disabled by default).<br /> 
- `--fastp` QC and quality assessment are performed with fastp instead of BBTools and FASTQC <br /> 
- `--genome` Set host genome. On the IKMB Medcluster valid options are human, mouse or chimp. In other cases this needs to be pre-configured. [How to add a host genome to the pipeline?](hostgenome.md) <br />
- `--no_qc` Skips QC-Module. Only use if your input reads are the output of `--cleanreads`<br /> 

## HUMAnN options
- `--metaphlan_db` Directory of Metaphlan database. REQUIRED! <br /> 
- `--humann_db`: Directory of HUMAnN database. REQUIRED! <br /> 

## Assembly options
- `--assemblymode` Specify assembly mode
    - **single** (default) Single-sample assembly
    - **group** Group-based co-assembly (requires input as CSV with `group` column).
    - **all** Cohort-wide co-assembly
> We recommend co-assembly with only moderate group sizes (~100 samples) due to hardware restrictions.<br />
- `--binner` Comma-separated list of binning tools (default: all). Options: **concoct**,**maxbin**,**semibin**,**metabat**,**vamb** <br />
- `--contigsminlength` Minimum contig length (default: 2000). <br />
- `--semibin_environment` Specify SemiBin2 environment (default: **human_gut**). See the [SemiBin Documentation](https://github.com/BigDataBiology/SemiBin/#easy-singleco-assembly-binning-mode) for other options. Choose **global** if no other environment is appropiate.  <br />
- `--skip_gtdbtk` Skip GTDB-TK for taxonomical assignment. <br />
- `--skip_checkm` Skip Checkm bin quality check. <br />
- `--gtdbtk_reference` Directory of GTDB-TK Reference.<br />
- `--publish_megahit` Publish assembled Megahit contigs.<br />
- `--publish_rawbins` Publish the results of all used binning tools in the genome assembly workflow.<br />
- `--vamb_groupsize` Set subgroup size for VAMB (default: 100). Adjust based on cohort size.<br />
### MAGScoT options
- `--magscot_min_sharing` Scoring parameter a [default=1] <br />
- `--magscot_score_a` Scoring parameter a [default=1] <br />
- `--magscot_score_b` Scoring parameter b [default=0.5] <br />
- `--magscot_score_c` Scoring parameter c [default=0.5] <br />
- `--magscot_threshold` Scoring minimum completeness threshold [default=0.5] <br />
- `--magscot_min_markers` Minimum number of unique markers in bins to be considered as seed for bin merging [default=25] <br />
- `--magscot_iterations` Number of merging iterations to perform. [default=2] <br />

## MetaPhlAn options
- `--metaphlan_db` Directory of Metaphlan database. REQUIRED! <br /> 
- `--publish_metaphlanbam` Publish the bam file output of Metaphlan. <br /> 

## Kraken2 options
- `--kraken2_db` Directory of used Kraken2 database. Should be Bracken ready for use with Bracken. REQUIRED! <br />

### Bracken options and their default
- `--bracken_length` = 100<br />
- `--bracken_level` = "S"<br />
- `--bracken_threshold` = 0<br />

## Sylph options
- `--sylph_db` Set the path to a sylph databse.<br />
- `--sylph_merge` All sylph profiling will be done in one process. Produces a single output for all samples combined. <br />
- `--sylph_processing` Shortcut for high-throughput data processing with sylph, skips quality control, no other modules available in this mode.  <br />

## Salmon options
>**Note**: The usage of Salmon for metagenomes is experimental.
- `--salmon_db` Directory of used salmon database. REQUIRED! <br />
- `--salmon_reference` Path to tab-separated taxonomy file corresponding to the used salmon database. Not required if used database contains taxonomic names. Two column file with header line containing in the first column the bin names used in the salmon database and in the second column the taxonomic assignment by GTDB-Tk in the format "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli". <br />



