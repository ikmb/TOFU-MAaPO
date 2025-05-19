#!/bin/bash
# Utility script to generate a CSV file which can be used as the input for the TOFU-MAaPO pipeline.
# Usage message
usage() {
    echo "Usage: $0 -i <input_dir> [-i <input_dir2> ...] -o <output.csv> [-r1 <suffix_read1>] [-r2 <suffix_read2>] [-u <suffix_unpaired>] [-e <file_ending>]"
    echo "Defaults: suffix_read1=_1, suffix_read2=_2, suffix_unpaired=_unpaired, file_ending=.fastq.gz"
    echo "Example: $0 -i dir1 -i dir2 -o output.csv -r1 _R1 -r2 _R2 -u _single -e .fq.gz"
    exit 1
}

# Default values
SUFFIX_R1="_1"
SUFFIX_R2="_2"
SUFFIX_UNPAIRED="_unpaired"
FILEENDING=".fastq.gz"
INPUT_DIRS=()

# Parse parameters
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            INPUT_DIRS+=("$(realpath "$2")")
            shift 2
            ;;
        -o|--output)
            OUTPUT="$(realpath "$2")"
            shift 2
            ;;
        -r1|--suffix-read1)
            SUFFIX_R1="$2"
            shift 2
            ;;
        -r2|--suffix-read2)
            SUFFIX_R2="$2"
            shift 2
            ;;
        -u|--suffix-unpaired)
            SUFFIX_UNPAIRED="$2"
            shift 2
            ;;
        -e|--file-ending)
            FILEENDING="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Ensure required args are provided
if [[ ${#INPUT_DIRS[@]} -eq 0 || -z "$OUTPUT" ]]; then
    usage
fi

# Append file ending to suffixes
SUFFIX_R1="$SUFFIX_R1$FILEENDING"
SUFFIX_R2="$SUFFIX_R2$FILEENDING"
SUFFIX_UNPAIRED="$SUFFIX_UNPAIRED$FILEENDING"

# Start CSV
echo "id,read1,read2,read3" > "$OUTPUT"

# Declare associative arrays
declare -A READ1 READ2 SINGLE UNPAIRED

# Collect files
ALL_FILES=()
for DIR in "${INPUT_DIRS[@]}"; do
    while IFS= read -r -d '' FILE; do
        ALL_FILES+=("$FILE")
    done < <(find "$DIR" -type f -name "*$FILEENDING" -print0)
done

# Sort collected files
IFS=$'\n' SORTED_FILES=($(sort <<<"${ALL_FILES[*]}"))
unset IFS

# Process files
for FILE in "${SORTED_FILES[@]}"; do
    [[ -e "$FILE" ]] || continue
    FILE_PATH=$(realpath "$FILE")
    BASENAME=$(basename "$FILE")

    if [[ $BASENAME =~ (.+)"$SUFFIX_R1"$ ]]; then
        ID="${BASH_REMATCH[1]}"
        READ1["$ID"]="$FILE_PATH"
    elif [[ $BASENAME =~ (.+)"$SUFFIX_R2"$ ]]; then
        ID="${BASH_REMATCH[1]}"
        READ2["$ID"]="$FILE_PATH"
    elif [[ $BASENAME =~ (.+)"$SUFFIX_UNPAIRED"$ ]]; then
        ID="${BASH_REMATCH[1]}"
        UNPAIRED["$ID"]="$FILE_PATH"
    elif [[ $BASENAME =~ (.+)\.fastq\.gz$ ]]; then
        ID="${BASH_REMATCH[1]}"
        if [[ -n "${READ1[$ID]}" || -n "${READ2[$ID]}" ]]; then
            UNPAIRED["$ID"]="$FILE_PATH"
        else
            SINGLE["$ID"]="$FILE_PATH"
        fi
    fi
done

# Write paired-end reads
for ID in "${!READ1[@]}"; do
    echo "$ID,${READ1[$ID]},${READ2[$ID]:-},${UNPAIRED[$ID]:-}" >> "$OUTPUT"
done

# Write single-end reads
for ID in "${!SINGLE[@]}"; do
    echo "$ID,${SINGLE[$ID]},," >> "$OUTPUT"
done

# Write orphan unpaired reads
for ID in "${!UNPAIRED[@]}"; do
    if [[ -z "${READ1[$ID]}" && -z "${READ2[$ID]}" ]]; then
        echo "$ID,,," >> "$OUTPUT"
    fi
done

echo "CSV file generated: $OUTPUT"
