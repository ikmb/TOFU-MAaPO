![](images/ikmb_bfx_logo.png)

# TOFU-MAaPO

Taxonomic Or FUnctional Metagenomic Assembly and PrOfiling = TOFU-MAaPO 

# Pipeline Structure
![](./images/metawo_overview.png)
Overview of TOFU-MAaPO 1.3.1

# Overview

The pipeline analyses metagenomic short reads and can perform analysis for taxonomic profiling, abundance of microbial metabolic pathways and assembly of metagenomic genomes. 


TOFU-MAaPO is a Nextflow pipeline for the processing and analysis of metagenomic short reads with a focus on gut metagenomes. The pipeline can run on any linux system and needs as dependencies only Nextflow and the container engine Singularity. An installation step is not needed as Nextflow can automatically download all needed files. <br />

The pipeline can install needed databases on its own, see for this the usage documentation. The pipeline can therefor be downloaded and executed in one single command line call. It requires short single- or paired-end metagenomic shotgun sequencing FASTQ files as input, or a single csv file containing a list of samples and their associated FASTQ files or can download data from SRA by giving the project or sample id as an input.<br />

TOFU-MAaPO processes the data for quality control and possible host decontamination (optional) and performs downstream analysis for taxonomic abundance profiles of each sample using Kraken2, Bracken, Sylph, Salmon and/or Metaphlan4, metabolic pathway analysis using HUMAnN (v3.6) and assembly of metagenomic genomes (MAGs).<br />

Genome assembly is done by generating contigs from the qc'ed reads with Megahit (single samples individually, grouped or all samples combined). The contigs are then catalogued and indexed using Minimap2 and then binned with the option to use up to five binning tools (Metabat2, Concoct, Maxbin, Semibin2 and vamb). The resulting bins will then be refined and, where possible, combined with MAGScoT based on sets of single-copy microbial marker genes from the Genome Taxonomy Database. The profiles of present marker genes in each result from the different binning algorithms are compared, and new hybrid candidate bins are created if the bin sets share a user-adjustable proportion of marker genes. The results are also taxonomically annotated with GTDB-TK and quality checked with Checkm. An estimated bin coverage per sample is generated as additional output. <br />

# Quick start
## Prerequisites:
This pipeline uses tools that require more computing capacity and memory than a workstation normally has. The tools Semibin for example can use up to 200GB of RAM, while GTDB-TK also uses 100GB of RAM to name two examples. Here, we show you how to run TOFU-MAaPO on a laptop with at least 4 cores and 32GB of memory, but we would not recommend this. A computer system for TOFU-MAaPO has at least 32 cores and a minimum of 128 GB of RAM. We recommend using it on an HPC for larger amounts of data. <br />

## Nextflow and Singularity installation
For Nextflow please install a suited Java version for Nextflow, we recommend to use SDKMAN for easy Java installation:
```
# Install Java Temurin with SDKMAN (other Java versions might cause bugs)
curl -s https://get.sdkman.io | bash
sdk install java 17.0.10-tem
#confirm that java is used in version 17.0.10-tem
java -version
#In case another java version is shown: Create and activate a sdk environment in the directory you want to execute the Nextflow pipeline
sdk env init
sdk env
```
Nextflow can be installed and tested with:
```
# Install Nextflow
curl -s https://get.nextflow.io | bash
# Make Nextflow executable:
chmod +x nextflow
# Try a simple Nextflow demo
nextflow run hello
```
Singularity (now Apptainer) can be install either by following the [Singularity Quickstart Guide](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html) or by installing it via
[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). With Conda, you can install Singularity without sudo rights. After installation, make sure that the environment is activated and test that Singularity is working:<br />
```
# Create a new conda environment for Singularity
conda create --name sing_env -c conda-forge -c bioconda singularity=3.8 
# Activate environment
conda activate sing_env
# Check whether Singularity has been successfully installed
singularity --version
# Also make sure you can run an example container
singularity run library://sylabsed/examples/lolcow
```

## Download or update TOFU-MAaPO
You can download and update the pipeline directly with Nextflow:<br />
```
nextflow pull ikmb/TOFU-MAaPO
```
You will find the pipeline code stored in `${HOME}/.nextflow/assets/ikmb/TOFU-MAaPO`.<br />


## Running TOFU-MAaPO:
The input for TOFU-MAaPO can be either locally stored fastq.gz files or SRA IDs (sample or project).<br />

### Locally stored input:
Create a directory to work in and download a publicly available paired-end metagenome  in your home-directory:
```
mkdir -p ${HOME}/tofu-quickstart && cd ${HOME}/tofu-quickstart
wget https://ibdmdb.org/downloads/raw/HMP2/MGX/2018-05-04/PSM6XBR1.tar
tar -xvf PSM6XBR1.tar && rm PSM6XBR1.tar
```

Now start TOFU-MAaPO to perform QC (without host read removal) on the example metagenome. This should only take some minutes. <br />
```
nextflow run ikmb/TOFU-MAaPO \
    -profile quickstart \
    --reads '*_R{1,2}.fastq.gz' \
    --cleanreads \
    --outdir results
```
With the parameter `--cleanreads` the qc'ed fastq files will be published into the output directory "results".

### SRA input:

TOFU-MAaPO can download metagenomes from NCBI SRA by SRA IDs. These can be Project, Sample or Run IDs. For this you need to get an account at NCBI and create your personal NCBI API key from NCBI. You can get your NCBI API key by going to NCBI -> Account -> [Account Settings](https://ncbi.nlm.nih.gov/account/settings/) -> API Key Management.

Create a directory to work in if you haven't done so earlier:
```
mkdir -p ${HOME}/tofu-quickstart && cd ${HOME}/tofu-quickstart
```

A workflow with QC (without host read removal) for a single metagenome with known SRA Run ID as input would look like this:

```
nextflow run ikmb/TOFU-MAaPO \
    -profile quickstart \
    ---sra 'SRX3105436' \
    --apikey **YOUR_NCBI_API_KEY** \
    --outdir results
```


### Okay, now we have run TOFU-MAaPO once for quality control, but I want to use some analysis tools of TOFU-MAaPO on my samples or want to run it on a large machine

Up to this point, we have used a profile/configuration of TOFU-MAaPO with the name "quickstart". This is only designed for TOFU-MAaPO to be run locally in the "home" directory, to use a maximum of 4 cores per process and to use only up to 32GB of memory.  
To use all the features of TOFU-MAaPO, it is necessary to edit a configuration file so that Nextflow knows, for example, whether we want to use a scheduler such as SLURM on an HPC or where the reference databases are or should be located. It is also possible to make far more memory and CPU cores available to the pipeline.
Please edit the configurations to your system needs as explained in the [installation and configuration documentation](docs/installation.md).<br />

On a computer system with more than 100GB of RAM and a configured custom.config file you can for example execute the examples from above but with additional MAG assembly. Beware, this might take some hours:
#### With locally available input:
```
nextflow run ikmb/TOFU-MAaPO \
    -profile custom \
    --reads '*_R{1,2}.fastq.gz' \
    --assembly \
    --updategtdbtk \
    --gtdbtk_reference '/path/to/download/gtdbtk_db/to' \
    --outdir results
```

#### Or with SRA input:
```
nextflow run ikmb/TOFU-MAaPO \
    -profile custom \
    ---sra 'SRX3105436' \
    --assembly \
    --updategtdbtk \
    --apikey **YOUR_NCBI_API_KEY**
    --gtdbtk_reference '/path/to/download/gtdbtk_db/to' \
    --outdir results
```

For further usage options please see the [usage documentation](docs/usage.md).<br />

## Configuration
Please edit the configurations to your system needs as explained in the [installation and configuration documentation](docs/installation.md).<br />

# Documentation 

All further documentation about the pipeline can be found in the `docs/` directory or under the links below:

1. [Overview/Flowchart of TOFU-MAaPO](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Available options](docs/usage.md)
4. [Output Structure](docs/output.md)


# Funding

The project was funded by the German Research Foundation (DFG) [Research Unit 5042 - miTarget INF](https://www.mitarget.org/).
