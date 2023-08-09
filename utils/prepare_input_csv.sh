#!/bin/bash

declare -A processed_files

while getopts "i:o:f:r:" opt; do
  case ${opt} in
    i ) input_dir=$OPTARG;;
    o ) output_file=$OPTARG;;
    f ) forward_suffix=$OPTARG;;
    r ) reverse_suffix=$OPTARG;;
    \? ) echo "Usage: cmd [-i input_dir] [-o output_file] [-f forward_suffix] [-r reverse_suffix]"
         exit 1;;
  esac
done

if [[ -z "$input_dir" || -z "$output_file" ]]; then
    echo "Usage: $0 -i <input_dir> -o <output_file> [-forward <forward_suffix>] [-reverse <reverse_suffix>]"
    exit 1
fi

# Set default suffixes if not provided
forward_suffix=${forward_suffix:-"_R1"}
reverse_suffix=${reverse_suffix:-"_R2"}

# Create CSV header
echo "id,read1,read2" > "$output_file"

# Find fastq files and write to CSV
for file in "$input_dir"/*.fastq*; do
    id=$(basename "$file" | cut -d'.' -f1 | sed -e "s/$forward_suffix.*//" | sed -e "s/$reverse_suffix.*//")

    read1=$(ls "$input_dir"/"$id"*"$forward_suffix"*.fastq.gz 2>/dev/null)
    read2=$(ls "$input_dir"/"$id"*"$reverse_suffix"*.fastq.gz 2>/dev/null)
# Make sure, only unique entries are omitted into the output file:
    if [[ -z "${processed_files[$id]}" ]]; then
      processed_files[$id]=1
      echo "$id,$read1,$read2" >> "$output_file"
    fi
done
