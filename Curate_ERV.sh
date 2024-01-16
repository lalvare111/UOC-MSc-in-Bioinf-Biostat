#!/bin/bash

################################################################################
# Curate_ERV.sh                                                                #
#                                                                              #
# The script systematically curates VCF files containing non-reference         # 
# TE insertions. It applies various filters, including sorting, indexing,      #
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
paths=("path/to/results") #example: ("WES/50x/TE_consensus" "WGS/50x/TE_consensus")

# path to ERVcaller folder
ERV_CALLER="/path/to/ERVcaller/"

# path to repeat masked bed file
REPEAT_MASKED_BED="/path/to/repeat/masked/file.bed/"

#path to bed file of the target exone sequencing regions (plus 200pb at both ends)
EXOME_BED="/path/to/exome/bed/file.bed" 
################################################################################

# Iterates in each parent directory
for path in "${paths[@]}"; do
    vcf_files=("$path"/*.vcf)
    for vcf_file in "${vcf_files[@]}"; do
    
    # Extract the filename without the path
    filename=$(basename "$vcf_file")

    # Remove the file extension
    filename_without_extension="${filename%.*}"
    
    # Filter ERVcaller files 
    bcftools sort "$vcf_file" -Oz -o "$vcf_file".gz
    bcftools index "$vcf_file".gz
    bcftools view -R "${EXOME_BED}" "$vcf_file".gz -o "${path}/${filename_without_extension}"_reg.vcf
    perl "${ERV_CALLER}"/Scripts/Filtering_VCF_v1.4.pl "${path}/${filename_without_extension}"_reg.vcf 9 0.2 3 | \
    bedtools intersect -v -a - -b "${REPEAT_MASKED_BED}" > "${path}/${filename_without_extension}"_filtered_nested.vcf
    echo "Found VCF file: $(basename "$vcf_file")"
    done
done
