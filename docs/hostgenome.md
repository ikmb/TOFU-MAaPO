# Host genome decontamination with TOFU-MAaPO - removal of host contamination reads

This tutorial explains how to remove host contamination reads from your data using the TOFU-MAaPO pipeline. To perform host decontamination, you need to configure the Bowtie2 index of your host genome in the pipeline settings. Below are two methods to achieve this.

---

## Option 1: Add a pre-indexed host genome to the pipeline (example: Human genome)

1. **Download Bowtie2 indexes**  
   	Download the Bowtie2 indexes for your host genome from sources like [AWS Indexes](https://benlangmead.github.io/aws-indexes/bowtie).
	>**Note**: Files must be unzipped and contain the `.bt2` suffix

2. **Configure the pipeline**  
   	Add the path to the basename of the Bowtie2 index files to your custom configuration file (e.g., `tofu.config`) before running the QC module.   
	Example configuration for the human genome:

	```groovy
	params {
		genomes {
			human {
				bowtie_index = "/path/to/your/references/iGenomes/references/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome"
			}
		}
	}
	```
3. **Run the pipeline**  
Use the custom profile (-profile custom) and specify the genome using the `--genome human` parameter to remove human reads from your data.
	```bash
	nextflow run ikmb/TOFU-MAaPO -profile custom -c tofu.config --reads '*_R{1,2}.fastq.gz' --genome human
	```

## Option 2: Add a host genome from NCBI (example: Wild boar genome)
To include a new host genome, download it in FASTA format and create Bowtie2 indexes as follows:
1. **Set up environment**  
	Create a new Conda environment with the necessary tools:

	```bash
	conda create --name=bowtie2 -c conda-forge -c bioconda bowtie2 ncbi-genome-download unzip

	conda activate bowtie2
	```

2. **Download the genome**  
	Search for the genome accession code on [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) and download the genome. Example: Wild boar genome (Sus scrofa, Sscrofa 11.1):
	```bash
	datasets download genome accession GCF_000003025.6 --include genome
	```
3. **Extract the genome**
	Unzip the downloaded file and rename the genome file:
	```bash
	unzip ncbi_dataset.zip
	# Change path to the new created directory containing the fasta:
	cd ncbi_dataset/data/GCF_000003025.6
	# rename the genome to genome.fna
	mv *.fna genome.fna
	```	
4. **Create Bowtie2 indexes**
   	```bash
	bowtie2-build genome.fna genome
	```	
> **Note**: Optionally, move the `genome.*` files to a dedicated directory for reference genomes.

5. **Configure the Pipeline**  
	Update your configuration file (e.g., tofu.config) with the following entry for the wild boar genome:
	```groovy
	params {
		genomes {
			boar {
				bowtie_index = "/path/to/the/Bowtie2Index/genome"
			}
		}
	}
	```

6. **Run the pipeline**
	```bash
	nextflow run ikmb/TOFU-MAaPO -profile custom -c tofu.config --reads '*_R{1,2}.fastq.gz' --genome boar
	```

---

For additional usage customization options, refer to the [TOFU-MAaPO usage documentation](usage.md)