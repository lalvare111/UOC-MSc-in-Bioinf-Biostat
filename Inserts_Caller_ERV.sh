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
    bam_name=$(find /home/ceromova/scratch_lui/MGE/samples -type f -name "${sample}*${technology}*${coverage}*.bam" -exec basename {} \; | sed 's/\.[^.]*$//')
    bam_path=$(find /home/ceromova/scratch_lui/MGE/samples -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \;)
    archivo_bam_result=$(find /home/ceromova/scratch_lui/MGE/samples -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \; | sed 's#/home/ceromova/scratch_lui/MGE/samples##')
    results_folder="/home/ceromova/scratch_lui/MGE/results2$archivo_bam_result/$name"
    
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
      perl /home/ceromova/scratch_lui/MGE/ERVcaller/ERVcaller_v1.4.pl \
        -i "$bam_name" \
        -f .bam \
        -H /home/ceromova/scratch_lui/MGE/ref/hs38DH.fa \
        -T /home/ceromova/scratch_lui/MGE/ERVcaller/Database/"$full_name" \
        -I "$bam_path"/ \
        -O "$results_folder/" \
        -t 10 -S 20 -BWA_MEM -G -d WES -r 150 >> "$log_file" 2>&1
      echo "Process end for ${sample}_${technology}_${coverage}: $(date)" >> "$log_file"
    fi
  done
done < "$1"  # Redirects the first input file to the 'while read' loop
