#!/bin/bash

# Check if the input files exist and are provided
if [ ! -f "$1" ] || [ ! -f "$2" ]; then
  echo "Input file(s) not found or provided."
  echo "Usage: $0 <input_file_1> <input_file_2>"
  exit 1
fi

# Read each line from the second input file containing names
mapfile -t names < "$2"

# Read each line from the first input file and process the information
while read -r line; do
  # Split each line into an array based on space delimiter
  read -ra values <<< "$line"
  
  # Check if the line contains at least three elements
  if [ "${#values[@]}" -lt 3 ]; then
    echo "Each line in the input file must contain at least three values: <sample> <technology> <coverage>"
    continue  # Skip to the next line if the format is incorrect
  fi
  
  # Extract values into variables
  sample="${values[0]}"
  technology="${values[1]}"
  coverage="${values[2]}"
  
  # Process each name from the second file for the current line in the first file
  for full_name in "${names[@]}"; do
    # Extract filename and extension
    filename=$(basename -- "$full_name")
    name="${filename%.*}"  # Extract filename without extension
    
    # Perform operations using the extracted values and current name
    echo "Sample: $sample, Technology: $technology, Coverage: $coverage, Database: $full_name"
  
    # Create necessary directories and perform processing based on the extracted values
    bam_name=$(find /media/luis/TOSHIBA_EXT/TFM/SAMPLES -type f -name "${sample}*${technology}*${coverage}*.bam" -exec basename {} \; | sed 's/\.[^.]*$//')
    bam_path=$(find /media/luis/TOSHIBA_EXT/TFM/SAMPLES -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \;)
    archivo_bam_result=$(find /media/luis/TOSHIBA_EXT/TFM/SAMPLES -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \; | sed 's#/media/luis/TOSHIBA_EXT/TFM/SAMPLES##')
    results_folder="/media/luis/TOSHIBA_EXT/TFM/Resultados$archivo_bam_result/$name"
    
    if [ -z "$bam_name" ]; then
      echo "No BAM file found for ${sample}_${technology}_${coverage}."
      continue  # Skip to the next name if BAM file not found
    else
      if [ ! -d "$results_folder" ]; then
        mkdir -p "$results_folder"
        echo "Directory created: $results_folder"
      else
        echo "Directory already exists: $results_folder"
      fi
      
      log_file="${name}_${sample}_${technology}_${coverage}.log"
      
      # Perform processing using extracted values (example: ERVcaller command)
      echo "Process start for ${sample}_${technology}_${coverage}: $(date)" > "$log_file"
      java -jar -Xmx4G /media/luis/TOSHIBA_EXT/TFM/MELTv2.2.2/MELT.jar Single \
      	-c $(echo "${coverage}" | sed 's/x//g') \
      	-h /media/luis/TOSHIBA_EXT/TFM/refs/hs38DH.fa \
      	-bamfile "${bam_path}/${bam_name}".bam \
      	-n /media/luis/TOSHIBA_EXT/TFM/MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed \
      	-t /media/luis/TOSHIBA_EXT/TFM/MELTv2.2.2/"$full_name" \
      	-w /media/luis/TOSHIBA_EXT/TFM/Resultados"$archivo_bam_result"/MELT_EX_TE_consensus/"${sample}" \
      	-exome true \
      	-e 150 >> "$log_file" 2>&1
      echo "Process end for ${sample}_${technology}_${coverage}: $(date)" >> "$log_file"
    fi
  done
done < "$1"  # Redirects the first input file to the 'while read' loop
