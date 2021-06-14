# Outputs

The pipeline generates the following outputs, located under "results", if no --outdir was set.

* A folder per sample with quality metrics and cleaned fastq files.
* `MultiQC/`A folder that contains the combined QC reports
* `pipeline_info` A folder containing workflow traces
* `humann` A folder that contains the merged result tables of HUMAnN3
* `Metaphlan3` A folder containing the merged abundancy table of Metaphlan3