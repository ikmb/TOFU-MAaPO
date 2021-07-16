# Usage:

This pipeline requires Nextflow 21.04.0 or higher. Other dependencies are containerized with Singularity and Docker.<br />

On Kiel Medcluster, please load the following modules with:
```bash
module load singularity nextflow
```

Reference databases for Metaphlan3, Kraken2 and HUMAnN3 are needed. On Kiel Medcluster, these are already set in the respective config file.
Metaphlan DB: `--metaphlan_db`<br />
HUMAnN DB:    `--humann_db`<br />
Kraken DB:    `--kraken2_db`<br />

Pipeline is module based and will run in the most basic run only QC steps.

Run the Pipeline with<br />
```bash
nextflow run ikmb/metagenomic-workflows --reads '/path/to/fastqfiles/*_R{1,2}_001.fastq.gz'
```
## Available modules:
For analysis following modules are available:<br />
**--virus** Run Kraken2 with a on Medcluster preconfigured RefSeq viruse database.<br />
**--metaphlan** Run Metaphlan3<br />
**--humann** Run HUMAnN3<br />

Experimental:<br />
**--assembly** Run a genome assembly workflow<br />
Reference database for GTDB-TK needs to be set (already set on Kiel Medcluster):<br />
GTDB-TK Reference: `--GTDBTKreference`<br />

## QC otions:
**--genome** set host genome. On the IKMB Medcluster valid options are human, mouse or chimp. In other cases this needs to be pre-configured.<br />

## Other options:
**--outdir** set a custom output directory, default is "results"<br />
**-resume** resumes pipeline and will continue the run with already completed, cached processes.<br />
