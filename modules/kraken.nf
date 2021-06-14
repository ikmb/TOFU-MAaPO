process KRAKEN2 {

label 'kraken'

publishDir "${params.outdir}/${sampleID}/Kraken/", mode: 'copy'

input:
tuple val(sampleID),file(left),file(right),file(unpaired)

output:
tuple val(sampleID),file(report), emit: krakenreport


script:
report = sampleID + ".kraken2_report.txt"
kraken_log = sampleID + ".kraken2.log"

"""
kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --output $kraken_log --report $report $left $right
 """
}
//output: tuple val(sampleID),file(kraken_log), emit: krakenlog

process KRAKEN2YAML {

input:
file(reports)

output:
file(report_yaml)

script:

report_yaml = "kraken_report_mqc.yaml"
"""	
kraken2yaml.pl --outfile $report_yaml
"""
}