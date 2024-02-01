# Outputs

The pipeline generates the following outputs depending on which modules were activated, located under "results", if no custom `--outdir` was set.

* `MultiQC`A folder that contains the combined QC reports
* `pipeline_info` A folder containing workflow traces and a file containing all software versions used in the performed workflow.
* `humann` A folder that contains the merged result tables of HUMAnN3
* `Metaphlan4` A folder containing the merged abundances tables of Metaphlan4
* `Kraken` A folder containing the merged abundances tables of Kraken
* `sylth` A folder containing the abundances profiles calculated with sylth
* `salmon` A folder containing the abundances profiles calculated with salmon
* `FastQC` A folder containing all FASTQC outputs as zip and html
Assembly outputs:
* `checkm` A folder containing all checkm outputs
* `GTDBTK` A folder containing all GTDB-TK outputs
* `MAG_abundance` A folder containing all estimated MAG abundances per sample
* `magscot` A folder containing all magscot outputs
* `counttable` A folder containing count tables for created contigs per sample