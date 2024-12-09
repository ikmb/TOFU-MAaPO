# Host genome decontamination with TOFU-MAaPO

For host decontamination you need to add the bowtie2 index of your host genome to the pipeline configurations. 
You have two options:
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