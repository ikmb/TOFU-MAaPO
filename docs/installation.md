# Installation

This pipeline is containerized with Singularity and Docker. You will need to prepare databases for respective modules (Metaphlan4, HUMAnN3 or Kraken2) and a config file for your compute system.

Download the pipeline into your home folder with:
```bash
nextflow pull ikmb/TOFU-MAaPO
```
The pipeline can then be found in `~/.nextflow/assets/ikmb/TOFU-MAaPO`

Edit the config file `conf/custom.config` to your system. Please follow the instructions below before running the respective modules within the pipeline. You will then be able to run the pipeline with your config with the parameter `-profile custom`. <br />

### Executors
Should you want to run the pipeline on a HPC or Cloud service, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html) to adapt your custom config file.

## QC

For host decontamination you need to add the bowtie2 index of your host genome to the pipeline configurations. Here are two solutions how you can do that:

### Add one already indexed host genome to the pipeline
Download your needed host genome as Bowtie2 indexes from e.g. [here](https://benlangmead.github.io/aws-indexes/bowtie) and set the the path to the  basename of the index files in your custom config file (we use the custom.config in this case) prior running the pipeline like so:
```
params {
	'genomes' {
		'human' { bowtie_index = "/path/to/your/references/iGenomes/references/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome"
		}
	}
}
```
You can now use with the custom profile `-profile custom` the parameter `--genome human` to remove human reads from your data.

### Add a host genome from NCBI to the pipeline

To add a host genome to your pipeline configurations, you need the genome in bowtie2-index format. For this, make sure, the genome is available in fasta format. In this example, we first download a genome from NCBI and create the index with tools installed with conda.
```bash
# Create a new conda environment for bowtie2 and the ncbi-genome-download tool
conda create --name=bowtie2 -c conda-forge -c bioconda bowtie2 ncbi-genome-download unzip
# Activate environment
conda activate bowtie2
# Search on ncbi.nlm.nih.gov/datasets for the accession code of the host genome of your choice, as an example we use pig (Sus scrofa) Sscrofa 11.1:
datasets download genome accession GCF_000003025.6 --include genome
# unzip the downloaded file
unzip ncbi_dataset.zip
# Change path to the new created directory containing the fasta:
cd ncbi_dataset/data/GCF_000003025.6
# rename the genome to genome.fna
mv *.fna genome.fna
# create bowtie2-index
bowtie2-build genome.fna genome
```

If you like, you can now move the genome.* files to a repository directory. Now edit your config file for TOFU-MAaPO (we use the custom.config in this case) to contain following entry:
```
params {
//reference genomes for host removal
	'genomes' {
		'pig' {
			bowtie_index = "/path/to/the/Bowtie2Index/genome"
		}
	}
}
```
In this example while using your edited config profile `-profile custom` you can now use `--genome pig` to remove pig read sequences from your data.

## Metaphlan & Humann
The pipeline needs databases for HUMAnN3, Metaphlan4 and Kraken2. For HUMAnN3 and Metaphlan4 the pipeline can download the necessary files with the parameters `--updatemetaphlan` and `--updatehumann` to the paths you supply with `--metaphlan_db` and `--humann_db`. 
Run the pipeline with the given parameters in your first run. The pipeline will then download and unzip the databases at this path. Make sure your local computer is connected with the internet. After the pipeline sucessfully finished the first run, change your custom config file to set metaphlan_db and humann_db to your now downloaded databases. Note that for running HUMAnN3 besides of the HUMAnN3 database you will also need the Metaphlan4 database (in Version vJan21).

Edit and include this snipped in your custom config:
```
params {
	metaphlan_db = "/path/to/your/databases/Metaphlan/4.0"
	humann_db = "/path/to/your/databases/Humann3/3.6"
}
```

## Kraken2
Download and extract a Kraken 2 database for example from [here](https://benlangmead.github.io/aws-indexes/k2) and add them to your custom config file or set a path each run with the parameter `--kraken2_db`. Make sure, this database is also Bracken ready if you want to use Bracken.

Edit and include this snipped in your custom config:
```
params {
	kraken2_db = "/path/to/your/databases/Kraken2/k2_viral_20210517"
}
```

## Genome assembly
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

No further action is required! <br />
Databases for all tools (with the exception of Salmon) and human, mouse and chimp as host genomes are already set.