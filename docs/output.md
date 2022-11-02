# Outputs

The pipeline generates the following outputs, located under "results", if no --outdir was set.

* A folder per sample with individual output files.
* `MultiQC`A folder that contains the combined QC reports
* `pipeline_info` A folder containing workflow traces
* `humann` A folder that contains the merged result tables of HUMAnN3
* `Metaphlan4` A folder containing the merged abundances tables of Metaphlan4
* `Kraken` A folder containing the merged abundances tables of Kraken
* `FastQC` A folder containing all FASTQC outputs as zip and html