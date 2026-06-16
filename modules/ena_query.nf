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
    set -euo pipefail

    ACCESSION_QUERY=\$(printf '%s' "${accession}" | tr -d '[:space:]')
    OUTPUT_FILE="ena_query_results.csv"
    API_URL="https://www.ebi.ac.uk/ena/portal/api/search"
    QUERY_CHUNK_SIZE=50

    # Split TSV manually so empty fields between tabs are preserved.
    split_tsv_line() {
        local line="\$1"
        local sep=\$'\\t'
        SPLIT_FIELDS=()

        while [[ "\$line" == *"\$sep"* ]]; do
            SPLIT_FIELDS+=("\${line%%\$sep*}")
            line="\${line#*\$sep}"
        done

        SPLIT_FIELDS+=("\$line")
    }

    field_value() {
        local index="\$1"

        if [[ "\$index" -ge 0 && "\$index" -lt "\${#SPLIT_FIELDS[@]}" ]]; then
            printf '%s' "\${SPLIT_FIELDS[\$index]}"
        fi
    }

    # Keep CSV output quote-free for downstream splitCsv; ENA values are not expected to contain commas.
    csv_escape() {
        local value="\$1"
        value="\${value//\\\"/}"
        printf '%s' "\$value"
    }

    # ENA separates multiple FASTQ URLs for a run with semicolons.
    split_fastq_ftp() {
        local line="\$1"
        local sep=";"
        FASTQ_ITEMS=()

        while [[ "\$line" == *"\$sep"* ]]; do
            FASTQ_ITEMS+=("\${line%%\$sep*}")
            line="\${line#*\$sep}"
        done

        FASTQ_ITEMS+=("\$line")
    }

    # Classify filenames using common paired/singleton labels, case-insensitively.
    fastq_item_type() {
        local item="\$1"
        local filename="\${item##*/}"
        local lower

        lower=\$(printf '%s' "\$filename" | tr '[:upper:]' '[:lower:]')

        case "\$lower" in
            *unpaired*|*single*|*singleton*|*orphan*)
                printf '%s' "single"
                return
                ;;
            *_1*|*-1*|*.1*|*_r1*|*-r1*|*.r1*|*_f1*|*-f1*|*.f1*|*_f.*|*-f.*|*.f.*|*_forward*|*-forward*|*.forward*|*_fwd*|*-fwd*|*.fwd*)
                printf '%s' "read1"
                return
                ;;
            *_2*|*-2*|*.2*|*_r2*|*-r2*|*.r2*|*_r.*|*-r.*|*.r.*|*_reverse*|*-reverse*|*.reverse*|*_rev*|*-rev*|*.rev*)
                printf '%s' "read2"
                return
                ;;
        esac

        printf '%s' "unknown"
    }

    # For triplets, ensure paired reads are first and singleton/unpaired reads are last.
    # If read1/read2 are detected but the third file is unlabeled, treat the remainder as singleton.
    reorder_fastq_ftp() {
        local fastq_ftp="\$1"
        local read_type
        local single_idx=-1
        local read1_idx=-1
        local read2_idx=-1
        local i
        local reordered=""

        if [[ -z "\$fastq_ftp" || "\$fastq_ftp" == "NA" ]]; then
            printf '%s' "\$fastq_ftp"
            return
        fi

        split_fastq_ftp "\$fastq_ftp"

        if [[ "\${#FASTQ_ITEMS[@]}" -ne 3 ]]; then
            printf '%s' "\$fastq_ftp"
            return
        fi

        for i in "\${!FASTQ_ITEMS[@]}"; do
            read_type=\$(fastq_item_type "\${FASTQ_ITEMS[\$i]}")

            case "\$read_type" in
                single)
                    if [[ "\$single_idx" -eq -1 ]]; then
                        single_idx="\$i"
                    fi
                    ;;
                read1)
                    if [[ "\$read1_idx" -eq -1 ]]; then
                        read1_idx="\$i"
                    fi
                    ;;
                read2)
                    if [[ "\$read2_idx" -eq -1 ]]; then
                        read2_idx="\$i"
                    fi
                    ;;
            esac
        done

        if [[ "\$single_idx" -eq -1 && "\$read1_idx" -ge 0 && "\$read2_idx" -ge 0 ]]; then
            for i in "\${!FASTQ_ITEMS[@]}"; do
                if [[ "\$i" -ne "\$read1_idx" && "\$i" -ne "\$read2_idx" ]]; then
                    single_idx="\$i"
                    break
                fi
            done
        fi

        if [[ "\$single_idx" -lt 0 || "\$single_idx" -eq 2 ]]; then
            printf '%s' "\$fastq_ftp"
            return
        fi

        for i in "\${!FASTQ_ITEMS[@]}"; do
            if [[ "\$i" -eq "\$single_idx" ]]; then
                continue
            fi

            if [[ -n "\$reordered" ]]; then
                reordered="\${reordered};"
            fi

            reordered="\${reordered}\${FASTQ_ITEMS[\$i]}"
        done

        printf '%s;%s' "\$reordered" "\${FASTQ_ITEMS[\$single_idx]}"
    }

    # Preserve one output row per missing accession if ENA returns no data for a chunk.
    append_missing_accessions() {
        local fallback_accessions="\$1"
        local missing_accession

        MISSING_ACCESSIONS=()
        IFS=',' read -r -a MISSING_ACCESSIONS <<< "\$fallback_accessions"

        for missing_accession in "\${MISSING_ACCESSIONS[@]}"; do
            if [[ -z "\$missing_accession" ]]; then
                continue
            fi

            printf '%s,NA,NA,NA,NA\\n' "\$(csv_escape "\$missing_accession")" >> "\$OUTPUT_FILE"
        done
    }

    # Parse one ENA TSV response by header name, then append rows in the fixed output order.
    append_response_rows() {
        local response_file="\$1"
        local fallback_accessions="\$2"
        local header
        local line
        local wrote_row=0

        if ! IFS= read -r header < "\$response_file"; then
            append_missing_accessions "\$fallback_accessions"
            return
        fi

        header="\${header%\$'\\r'}"
        split_tsv_line "\$header"

        # Initialize indexes so missing columns cleanly become NA instead of shifting values.
        RUN_ACCESSION_IDX=-1
        FASTQ_FTP_IDX=-1
        FIRST_PUBLIC_IDX=-1
        INSTRUMENT_MODEL_IDX=-1
        READ_COUNT_IDX=-1

        for i in "\${!SPLIT_FIELDS[@]}"; do
            case "\${SPLIT_FIELDS[\$i]}" in
                run_accession) RUN_ACCESSION_IDX="\$i" ;;
                fastq_ftp) FASTQ_FTP_IDX="\$i" ;;
                first_public) FIRST_PUBLIC_IDX="\$i" ;;
                instrument_model) INSTRUMENT_MODEL_IDX="\$i" ;;
                read_count) READ_COUNT_IDX="\$i" ;;
            esac
        done

        while IFS= read -r line || [[ -n "\$line" ]]; do
            line="\${line%\$'\\r'}"

            if [[ -z "\$line" ]]; then
                continue
            fi

            split_tsv_line "\$line"

            RUN_ACCESSION="\$(field_value "\$RUN_ACCESSION_IDX")"
            FASTQ_FTP="\$(field_value "\$FASTQ_FTP_IDX")"
            FIRST_PUBLIC="\$(field_value "\$FIRST_PUBLIC_IDX")"
            INSTRUMENT_MODEL="\$(field_value "\$INSTRUMENT_MODEL_IDX")"
            READ_COUNT="\$(field_value "\$READ_COUNT_IDX")"

            # Empty ENA fields are explicit NA values in the CSV.
            [[ -z "\$RUN_ACCESSION" ]] && RUN_ACCESSION="NA"
            [[ -z "\$FASTQ_FTP" ]] && FASTQ_FTP="NA"
            [[ -z "\$FIRST_PUBLIC" ]] && FIRST_PUBLIC="NA"
            [[ -z "\$INSTRUMENT_MODEL" ]] && INSTRUMENT_MODEL="NA"
            [[ -z "\$READ_COUNT" ]] && READ_COUNT="NA"
            FASTQ_FTP="\$(reorder_fastq_ftp "\$FASTQ_FTP")"

            printf '%s,%s,%s,%s,%s\\n' \\
                "\$(csv_escape "\$RUN_ACCESSION")" \\
                "\$(csv_escape "\$FASTQ_FTP")" \\
                "\$(csv_escape "\$FIRST_PUBLIC")" \\
                "\$(csv_escape "\$INSTRUMENT_MODEL")" \\
                "\$(csv_escape "\$READ_COUNT")" >> "\$OUTPUT_FILE"

            wrote_row=1
        done < <(tail -n +2 "\$response_file")

        if [[ "\$wrote_row" -eq 0 ]]; then
            append_missing_accessions "\$fallback_accessions"
        fi
    }

    # Normalize and validate the comma-separated accession list supplied by Nextflow.
    if [[ -z "\$ACCESSION_QUERY" ]]; then
        echo "Error: No accession provided!"
        exit 1
    fi

    IFS=',' read -r -a RAW_ACCESSIONS <<< "\$ACCESSION_QUERY"
    ACCESSIONS=()

    for QUERY_ACCESSION in "\${RAW_ACCESSIONS[@]}"; do
        if [[ -z "\$QUERY_ACCESSION" ]]; then
            continue
        fi

        ACCESSIONS+=("\$QUERY_ACCESSION")
    done

    if [[ "\${#ACCESSIONS[@]}" -eq 0 ]]; then
        echo "Error: No valid accession provided!"
        exit 1
    fi

    echo "run_accession,fastq_ftp,first_public,instrument_model,read_count" > "\$OUTPUT_FILE"

    # ENA queries can become too large, so send accessions in fixed-size chunks
    # and merge all response rows into the single output CSV.
    CHUNK_START=0
    CHUNK_ID=0

    while [[ "\$CHUNK_START" -lt "\${#ACCESSIONS[@]}" ]]; do
        CHUNK_END=\$((CHUNK_START + QUERY_CHUNK_SIZE))

        if [[ "\$CHUNK_END" -gt "\${#ACCESSIONS[@]}" ]]; then
            CHUNK_END="\${#ACCESSIONS[@]}"
        fi

        ENA_QUERY=""
        CHUNK_ACCESSIONS=""

        # Build one ENA search expression per chunk:
        # run_accession="SRR1" OR run_accession="SRR2" ...
        for ((i = CHUNK_START; i < CHUNK_END; i++)); do
            QUERY_ACCESSION="\${ACCESSIONS[\$i]}"

            if [[ -n "\$ENA_QUERY" ]]; then
                ENA_QUERY="\${ENA_QUERY} OR "
                CHUNK_ACCESSIONS="\${CHUNK_ACCESSIONS},"
            fi

            ENA_QUERY="\${ENA_QUERY}run_accession=\\\"\${QUERY_ACCESSION}\\\""
            CHUNK_ACCESSIONS="\${CHUNK_ACCESSIONS}\${QUERY_ACCESSION}"
        done

        RESPONSE_FILE="ena_response_\${CHUNK_ID}.tsv"
        curl -G -fsSL --retry 3 \\
            --data-urlencode "result=read_run" \\
            --data-urlencode "fields=run_accession,read_count,fastq_ftp,instrument_model,first_public" \\
            --data-urlencode "format=tsv" \\
            --data-urlencode "query=\$ENA_QUERY" \\
            "\$API_URL" > "\$RESPONSE_FILE"

        append_response_rows "\$RESPONSE_FILE" "\$CHUNK_ACCESSIONS"

        CHUNK_START="\$CHUNK_END"
        CHUNK_ID=\$((CHUNK_ID + 1))
    done

    """
}
