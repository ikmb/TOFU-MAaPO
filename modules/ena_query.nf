process ena_query {
	label 'local_download'
	scratch params.scratch
	publishDir "${params.outdir}/metadata", mode: 'copy', pattern: "*.csv"

	input:
	val(accession)

	output:
	path("ena_query_results.csv"), optional: true

	script:    
    """
    #!/bin/bash
    echo "${accession}" | awk -F, '{for(i=1;i<=NF;i++) print \$i}' > accession_input.csv

    INPUT_FILE="accession_input.csv"
    OUTPUT_FILE="ena_query_results_unfiltered.csv"
    API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&fields=read_count,fastq_ftp,instrument_model,first_public&accession=" 

    # Check if input file exists
    if [[ ! -f "\$INPUT_FILE" ]]; then
        echo "Error: Input file '\$INPUT_FILE' not found!"
        exit 1
    fi

    # Create output file with header if it doesn't exist
    if [[ ! -f "\$OUTPUT_FILE" ]]; then
        echo "run_accession,fastq_ftp,first_public,instrument_model,read_count" > "\$OUTPUT_FILE"
    fi

    # Process each ID in input file
    while IFS= read -r ID; do
        if [[ -z "\$ID" ]]; then
            continue  # Skip empty lines
        fi

        # Query API
        RESPONSE=\$(curl -s "\${API_URL}\${ID}")
        # Skip the header row from response and extract values
        echo "\$RESPONSE" | tail -n +2 | while IFS=\$'\t' read -r RUN_ACCESSION FASTQ_FTP FIRST_PUBLIC INSTRUMENT_MODEL READ_COUNT; do
            # Handle empty values by replacing with "N/A"
            [[ -z "\$RUN_ACCESSION" ]] && RUN_ACCESSION="NA"
            [[ -z "\$FASTQ_FTP" ]] && FASTQ_FTP="NA"
            [[ -z "\$FIRST_PUBLIC" ]] && FIRST_PUBLIC="NA"
            [[ -z "\$INSTRUMENT_MODEL" ]] && INSTRUMENT_MODEL="NA"
            [[ -z "\$READ_COUNT" ]] && READ_COUNT="NA"

            # Append the result to CSV
            echo "\$RUN_ACCESSION,\$FASTQ_FTP,\$FIRST_PUBLIC,\$INSTRUMENT_MODEL,\$READ_COUNT" >> "\$OUTPUT_FILE"
        done
    done < "\$INPUT_FILE"

    # Remove entries that are just looking like the header
    awk 'NR == 1 || \$0 !~ /run_accession,fastq_ftp,first_public,instrument_model,read_count/' \$OUTPUT_FILE > ena_query_results.csv

    """
}