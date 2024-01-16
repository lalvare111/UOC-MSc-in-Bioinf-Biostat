#!/bin/bash

################################################################################
# Curate_MELT.sh                                                               #
#                                                                              #
# The script systematically curates VCF files containing non-reference         # 
# TE insertions. Concatenates the different MELT output files and              # 
# applies various filters, including sorting, indexing,                        #
# and filtering based on read annotations and specified BED files.             #   
#                                                                              #
# author: Luis Alvarez Fernandez                                               #
# version: 1.0                                                                 #
# email: lalvarezfernandez@uoc.edu                                             #
# academic year: 2023-2014                                                     #
################################################################################

###############################################################################
# CONFIGURATION
# Defines an array of the main directories
paths=("path/to/results") #example: (WES/50x/MELT_EX_TE_consensus" "WGS/50x/MELT_TE_consensus")

# path to repeat masked bed file
REPEAT_MASKED_BED="/path/to/repeat/masked/file.bed/"

#path to bed file of the target exone sequencing regions (plus 200pb at both ends)
EXOME_BED="/path/to/exome/bed/file.bed" 
################################################################################

# Iterates in each parent directory
for path in "${paths[@]}"; do
    # Finds and applies tabix -p vcf to each .vcf file
    find "$path" -type f -name "*.vcf" -exec sh -c '
        for vcf_file; do
            # Check VCF file size
            file_size=$(wc -c < "$vcf_file")
            if [ "$file_size" -gt 0 ]; then
                # Compress the VCF file with bgzip and move it to the output directory.
                bgzip -c "$vcf_file" > "${vcf_file}.gz"

                # Create an index file using tabix
                tabix -p vcf "${vcf_file}.gz"
            fi
        done' sh {} +

    # Find the subdirectories and concatenate the generated files.
    for sub_dir in "$path"/*; do
        # Check if it is a directory
        if [ -d "$sub_dir" ]; then
            # Concatenated file name
            concatenated_file_name=$(basename "$sub_dir")_"$(basename "$(dirname "$path")")_concat"

            # Concatenates .vcf.gz files in each subdirectory
            # Filter positions
            bcftools concat -a "$sub_dir"/*.vcf.gz -o "$sub_dir/$concatenated_file_name".vcf.gz
            bcftools index "$sub_dir/$concatenated_file_name".vcf.gz
            bcftools view -R "${EXOME_BED}" "$sub_dir/$concatenated_file_name".vcf.gz -o "$sub_dir/$concatenated_file_name"_reg.vcf
            bcftools view -i 'GT!="0/0" && DP>=6 && AD>=2 && (AD/(DP-AD))>=0.2' "$sub_dir/$concatenated_file_name"_reg.vcf | \
            bedtools intersect -v -a - -b "${REPEAT_MASKED_BED}" > "$sub_dir/$concatenated_file_name"_filtered_nested.vcf
        fi
    done
done
