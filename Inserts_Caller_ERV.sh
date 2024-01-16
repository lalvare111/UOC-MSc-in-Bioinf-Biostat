#!/bin/bash

################################################################################
# Inserts_Caller_ERV.sh                                                        #
#                                                                              #
# The script finds the BAM files for the selected samples, executes            #
# the Perl script ERVcaller_v1.4.pl with the specified arguments, stores       #
# the output in the appropriate directories and keeps a log file               #
# of the process.                                                              #
#                                                                              #
# author: Luis Alvarez Fernandez                                               #
# version: 1.0                                                                 #
# email: lalvarezfernandez@uoc.edu                                             #
# academic year: 2023-2014                                                     #
################################################################################

###############################################################################
# CONFIGURATION
# path to ERVcaller folder
ERV_CALLER="/path/to/ERVcaller/"

#path to samples folder
SAMPLES_FOLDER="/path/to/samples/"

#path to reference genome
REF_GENOME="/path/to/ref/"

#ERV parameters
#Specify the number of threads
THREADS="number"

#Specify the length of split reads used (20 bp by default; >=40 bp is recommend for reads in length of 150 bp)
SPLIT_LENGHT="number"

#Specify read length (bp), including: 100, 150, and 250 bp
READ_LEN="number"
################################################################################
 
# Check if the input files exist and are provided
if [ ! -f "$1" ] || [ ! -f "$2" ]; then
  echo "Input file(s) not found or provided."
  echo "Usage: $0 <samples_file> <database_file>"
  exit 1
fi

# Read each line from the database_file containing database names
mapfile -t names < "$2"

# Read each line from the samples_file and process the information
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
  
  # Process each name from the database_file for the current line in the samples_file
  for full_name in "${names[@]}"; do
    # Extract filename and extension
    filename=$(basename -- "$full_name")
    name="${filename%.*}"  # Extract filename without extension
    
    # Perform operations using the extracted values and current name
    echo "Sample: $sample, Technology: $technology, Coverage: $coverage, Database: $full_name"
  
    # Create necessary directories and perform processing based on the extracted values
    bam_name=$(find "${SAMPLES_FOLDER}" -type f -name "${sample}*${technology}*${coverage}*.bam" -exec basename {} \; | sed 's/\.[^.]*$//')
    bam_path=$(find "${SAMPLES_FOLDER}" -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \;)
    archivo_bam_result=$(find "${SAMPLES_FOLDER}" -type f -name "${sample}*${technology}*${coverage}*.bam" -exec dirname {} \; | sed 's#${SAMPLES_FOLDER}##')
    results_folder="${SAMPLES_FOLDER}/results$archivo_bam_result/$name"
    
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
      
      # Perform processing using extracted values: ERVcaller script
      echo "Process start for ${sample}_${technology}_${coverage}: $(date)" > "$log_file"
      perl "${ERV_CALLER}"/ERVcaller_v1.4.pl \
        -i "$bam_name" \
        -f .bam \
        -H "${REF_GENOME}" \
        -T "${ERV_CALLER}"/Database/"$full_name" \
        -I "$bam_path"/ \
        -O "$results_folder/" \
        -t "${THREADS}" -S "${SPLIT_LENGHT}" -BWA_MEM -G -r "${READ_LEN}" >> "$log_file" 2>&1
      echo "Process end for ${sample}_${technology}_${coverage}: $(date)" >> "$log_file"
    fi
  done
done < "$1"  # Redirects the first input file to the 'while read' loop
