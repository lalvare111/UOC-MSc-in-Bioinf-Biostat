############################CASE SAMPLES: DIABETES (D)########################################################################
# File load
count_alleles_d <- read.table("Count_alleles_D")

# Create a new variable with the insertion type
count_alleles_d$Class <- gsub("<INS_MEI:(LINE1|ALU|SVA|HERV)>", "\\1", count_alleles_d$V5)

# Rename the variables of interest
colnames(count_alleles_d)[colnames(count_alleles_d) == "V1"] <- "chr"
colnames(count_alleles_d)[colnames(count_alleles_d) == "V2"] <- "pos"
colnames(count_alleles_d)[colnames(count_alleles_d) == "V10"] <- "het"
colnames(count_alleles_d)[colnames(count_alleles_d) == "V11"] <- "hom"

# Calculate relative and absolute allele frequencies
count_alleles_d$freq <- ((count_alleles_d$het + 2*count_alleles_d$hom)*100)/298
count_alleles_d$freq2 <- count_alleles_d$het + 2*count_alleles_d$hom

# Represent the allele frequencies
library(ggplot2)
ggplot(count_alleles_d, aes(x=freq, fill=Class)) + 
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.8, color = "black") +
  facet_wrap(~Class, scales = "fixed", ncol = 1) +
  labs(title = "Allele distribution by insertion type", x = "Frecuency", y = "Number of observations") +
  theme_classic() +
  theme(legend.position = "none")
  
# Represent the distribution by chromosomes
# Specify the order of the chromosomes and transform it into a factor
count_alleles_d$chr <- factor(count_alleles_d$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

ggplot(count_alleles_d, aes(x = chr, y = freq, color = chr)) +
  geom_point(size = 3) +
  labs(title = "Chromosomal distribution of insertion alleles", x = "Chromosome", y = "Frequency") +
  theme_classic() +
  theme(legend.position = "none")

############################CONTROL SAMPLES: NON-DIABETES (ND)########################################################################
# File load
count_alleles_nd <- read.table("Count_alleles_ND")

# Create a new variable with the insertion type
count_alleles_nd$Class <- gsub("<INS_MEI:(LINE1|ALU|SVA|HERV)>", "\\1", count_alleles_nd$V5)

# Rename the variables of interest
colnames(count_alleles_nd)[colnames(count_alleles_nd) == "V1"] <- "chr"
colnames(count_alleles_nd)[colnames(count_alleles_nd) == "V2"] <- "pos"
colnames(count_alleles_nd)[colnames(count_alleles_nd) == "V10"] <- "het"
colnames(count_alleles_nd)[colnames(count_alleles_nd) == "V11"] <- "hom"

# Calculate relative and absolute allele frequencies
count_alleles_nd$freq <- ((count_alleles_nd$het + 2*count_alleles_nd$hom)*100)/1178
count_alleles_nd$freq2 <- count_alleles_nd$het + 2*count_alleles_nd$hom

# Represent the allele frequencies
library(ggplot2)
ggplot(count_alleles_nd, aes(x=freq, fill=Class)) + 
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.8, color = "black") +
  facet_wrap(~Class, scales = "fixed", ncol = 1) +
  labs(title = "Allele distribution by insertion type", x = "Frecuency", y = "Number of observations") +
  theme_classic() +
  theme(legend.position = "none")

# Represent the distribution by chromosomes
# Specify the order of the chromosomes and transform it into a factor
count_alleles_nd$chr <- factor(count_alleles_nd$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

ggplot(count_alleles_nd, aes(x = chr, y = freq, color = chr)) +
  geom_point(size = 3) +
  labs(title = "Chromosomal distribution of insertion alleles", x = "Chromosome", y = "Frequency") +
  theme_classic() +
  theme(legend.position = "none")

##########################GENE ANOTATION###################################################################
BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)

# Create the folder if it doesn't exist to save the annotation file
system("[ ! -d 'BED_files' ] && mkdir 'BED_files'")

# Download the gene gtf file corresponding to the latest version () and transform it into a bed file
system("wget -O - ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz | gunzip -c - > BED_files/gencodeGenes.gtf")

# Import the annotation file
gencode_annotation <- makeTxDbFromGFF("BED_files/gencodeGenes.gtf", format = "gtf")


#################################CASE SAMPLES (D)#########################################################################3
# Filter insertions with occurrences in just one sample
count_alleles_d_filter <- count_alleles_d %>% filter(!(het %in% c(0) & hom %in% c(1)) & !(het %in% c(1) & hom %in% c(0)))

# Create a GRanges object from your dataonly
granges_d <- GRanges(seq = count_alleles_d_filter$chr, IRanges(start = count_alleles_d_filter$pos, width = 1))

# Find overlaps with gene_annotation
overlaps_d <- findOverlaps(granges_d, genes(gencode_annotation))

# Create a data frame with the hits 
overlaps_d_df <- data.frame(queryHits = queryHits(overlaps_d), subjectHits = subjectHits(overlaps_d))

# We keep just one hit per query (firts appearence)
overlaps_d_uniq <- overlaps_d_df %>%
  group_by(queryHits) %>%
  slice(1) %>%
  ungroup()

# Recover the gene annotation
annotation_data_d <- data.frame(ENS_gene=mcols(genes(gencode_annotation))[overlaps_d_uniq$subjectHits, ])

# Incorporate the annotation to the filtered data
count_alleles_d_filter$ENS_gene <- annotation_data_d$ENS_gene

# Add to the annotation of the genes the common name
library(biomaRt)
library(dbplyr)
library(BiocFileCache)

# Load the corresponding database
ensembl_dataset <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Save the external gene name in a new variable
attributes_vector_d <- character(length(count_alleles_d_filter$ENS_gene))

for (i in seq_along(attributes_vector_d)) {
  # Check if the value is already retrieved
  if (attributes_vector_d[i] == "") {
    # Get external_gene_name for each value
    gene_info <- getBM(
      attributes = "external_gene_name",
      filters = "ensembl_gene_id_version",
      values = count_alleles_d_filter$ENS_gene[i],
      mart = ensembl_dataset
    )
    # Assign external_gene_name to the corresponding position in the vector
    attributes_vector_d[i] <- gene_info$external_gene_name
  }
} 

count_alleles_d_filter$external_gene_name <- attributes_vector_d

#################################CONTROL SAMPLES (ND)#########################################################################3
# Filter insertions with occurrences in just one sample
count_alleles_nd_filter <- count_alleles_nd %>% filter(!(het %in% c(0) & hom %in% c(1)) & !(het %in% c(1) & hom %in% c(0)))

# Create a GRanges object
granges_nd <- GRanges(seq = count_alleles_nd_filter$chr, IRanges(start = count_alleles_nd_filter$pos, width = 1))

# Find overlaps with gene_annotation
overlaps_nd <- findOverlaps(granges_nd, genes(gencode_annotation))

# Create a data frame with the hits 
overlaps_nd_df <- data.frame(queryHits = queryHits(overlaps_nd), subjectHits = subjectHits(overlaps_nd))

# We keep just one hit per query (firts appearence)
overlaps_nd_uniq <- overlaps_nd_df %>%
  group_by(queryHits) %>%
  slice(1) %>%
  ungroup()

# Recover the gene annotation
annotation_data_nd <- data.frame(ENS_gene=mcols(genes(gencode_annotation))[overlaps_nd_uniq$subjectHits, ])

# We incorporate the annotation to the filtered data
count_alleles_nd_filter$ENS_gene <- annotation_data_nd$ENS_gene

# We add to the annotation of the genes the common name
attributes_vector_nd <- character(length(count_alleles_nd_filter$ENS_gene))

for (i in seq_along(attributes_vector_nd)) {
  # Check if the value is already retrieved
  if (attributes_vector_nd[i] == "") {
  # Get external_gene_name for each value
  gene_info <- getBM(
    attributes = "external_gene_name",
    filters = "ensembl_gene_id_version",
    values = count_alleles_nd_filter$ENS_gene[i],
    mart = ensembl_dataset
  )
  # Assign external_gene_name to the corresponding position in the vector
  attributes_vector_nd[i] <- gene_info$external_gene_name
  }
} 

count_alleles_nd_filter$external_gene_name <- attributes_vector_nd

##############################MERGE DATABASES#######################################################
# Keep the fields of interest in the databases
T2D_samples <-data.frame(ID = paste(count_alleles_d_filter$chr, count_alleles_d_filter$pos, sep=":"), count_alleles_d_filter[, c("Class", "freq", "hom", "het", "external_gene_name")])
CONTROL_samples <-data.frame(ID = paste(count_alleles_nd_filter$chr, count_alleles_nd_filter$pos, sep=":"), count_alleles_nd_filter[, c("Class", "freq", "hom", "het", "external_gene_name")])

Merge <- merge(T2D_samples, CONTROL_samples, by.x = 1, by.y = 1, all = TRUE)

library(gtools)
Merge_sort <- Merge[mixedorder(Merge$ID), ]

# Sort the data and convert the number of heterozygotes and homozygotes into frequencies
Merge_sort$freq.x <- round(Merge_sort$freq.x,1)
Merge_sort$freq.y <- round(Merge_sort$freq.y,1)
Merge_sort$hom.x <- round((Merge_sort$hom.x/149*100),1)
Merge_sort$hom.y <- round((Merge_sort$hom.y/598*100),1)
Merge_sort$het.x <- round((Merge_sort$het.x/149*100),1)
Merge_sort$het.y <- round((Merge_sort$het.y/598*100),1)

# Replace the NAs by 0s
Merge_sort <- Merge_sort %>%
  mutate_all(~ifelse(is.na(.), 0, .))

# Create a new variable with all the detected genes
Merge_sort <- mutate(Merge_sort, external_gene_name = ifelse(external_gene_name.x == 0, external_gene_name.y, external_gene_name.x))

# Merge those variables that share the same gene
Merge_sort2 <- Merge_sort %>%
  group_by(external_gene_name) %>%
  summarize(across(c(ID, Class.x, freq.x, hom.x, het.x, Class.y, freq.y, hom.y, het.y),
                   ~ ifelse(.x[1] != 0, as.character(.x[1]), as.character(.x[2])),
                   .names = "{.col}"),
            .groups = 'drop')

# Replace NAs by 0s
Merge_sort2 <- Merge_sort2 %>%
  mutate_all(~ifelse(is.na(.), 0, .))

# Stack the data in order to be able to represent them
Merge_sort_TD2 <- Merge_sort2[ ,c("ID", "freq.x", "het.x", "hom.x", "external_gene_name")]
colnames(Merge_sort_TD2) <- c("ID","freq","het","hom","gene_name")
Merge_sort_TD2$Sample_type <- rep("case", nrow(Merge_sort_TD2))
Merge_sort_TD2$freq <- as.numeric(Merge_sort_TD2$freq)

Merge_sort_Control <- Merge_sort2[ ,c("ID", "freq.y", "het.y", "hom.y", "external_gene_name")]
colnames(Merge_sort_Control) <- c("ID","freq","het","hom","gene_name")
Merge_sort_Control$Sample_type <- rep("control", nrow(Merge_sort_Control))
Merge_sort_Control$freq <- - as.numeric(Merge_sort_Control$freq)

Merge <- rbind(Merge_sort_TD2,Merge_sort_Control)

# Re-order the data frame according to chr:pos 
Merge_sort3 <- Merge[mixedorder(as.character(Merge$ID)), ]

# We define the variable ID and gene_name as factors (specifying the levels) as we will use them as labels on the y-axis scale
Merge_sort3$gene_name <- factor(Merge_sort3$gene_name, levels = rev(unique(Merge_sort3$gene_name)))
Merge_sort3$ID <- factor(Merge_sort3$ID, levels = rev(unique(Merge_sort3$ID)))

# We define the scale of the x-axis
freq_range <- range(as.numeric(Merge_sort3$freq))
freq_range_seq <- seq(freq_range[1], freq_range[2], by = 10)

library(ggh4x)

# Create the population pyramid type column chart for the allele frequencies of the insetions
g <- ggplot(Merge_sort3,
       aes(x = freq,
           y = gene_name, fill=Sample_type)) + 
  geom_col() + 
  theme_classic() + 
  labs(title = "Allele frequency of insertion per gene in the case and control samples", x = "Frequency", y = "Gene name") +
  scale_x_continuous(breaks  = freq_range_seq,
                     labels = abs(round(freq_range_seq,0))) 

# Add a second scale to the y-axis to display the chr:pos information
g + guides(y.sec = guide_axis_manual(
  labels = levels(Merge_sort3$ID)))

# create the population pyramid type column chart for homozygote and heterozygote frequencies
# Separate hom and het data into cases and controls to give them form

Merge_sort_TD2_het <- Merge_sort_TD2[, c("ID", "het", "gene_name", "Sample_type")]
Merge_sort_TD2_het$Sample_type <- rep("case_het", nrow(Merge_sort_TD2_het))
colnames(Merge_sort_TD2_het) <- c("ID", "freq", "gene_name", "Sample_type")
Merge_sort_TD2_het$freq <- as.numeric(Merge_sort_TD2_het$freq)

Merge_sort_TD2_hom <- Merge_sort_TD2[, c("ID", "hom", "gene_name", "Sample_type")]
Merge_sort_TD2_hom$Sample_type <- rep("case_hom", nrow(Merge_sort_TD2_hom))
colnames(Merge_sort_TD2_hom) <- c("ID", "freq", "gene_name", "Sample_type")
Merge_sort_TD2_hom$freq <- as.numeric(Merge_sort_TD2_hom$freq)

Merge_sort_Control_het <- Merge_sort_Control[, c("ID", "het", "gene_name", "Sample_type")]
Merge_sort_Control_het$Sample_type <- rep("control_het", nrow(Merge_sort_Control_het))
colnames(Merge_sort_Control_het) <- c("ID", "freq", "gene_name", "Sample_type")
Merge_sort_Control_het$freq <- - as.numeric(Merge_sort_Control_het$freq)

Merge_sort_Control_hom <- Merge_sort_Control[, c("ID", "hom", "gene_name", "Sample_type")]
Merge_sort_Control_hom$Sample_type <- rep("control_hom", nrow(Merge_sort_Control_hom))
colnames(Merge_sort_Control_hom) <- c("ID", "freq", "gene_name", "Sample_type")
Merge_sort_Control_hom$freq <- - as.numeric(Merge_sort_Control_hom$freq)

Merge2 <- rbind(Merge_sort_TD2_het,Merge_sort_TD2_hom,Merge_sort_Control_het,Merge_sort_Control_hom)

# Re-order the data frame according to chr:pos 
Merge_sort4 <- Merge2[mixedorder(as.character(Merge2$ID)), ]

# define the variable ID and gene_name as factors (specifying the levels) as we will use them as labels on the y-axis scale
Merge_sort4$gene_name <- factor(Merge_sort4$gene_name, levels = rev(unique(Merge_sort4$gene_name)))
Merge_sort4$ID <- factor(Merge_sort4$ID, levels = rev(unique(Merge_sort4$ID)))

library(ggh4x)
# Define the customised colours
custom_colors <- c("#FF9999", "#cc3333", "#66ccff", "#006699")


# Create the population pyramid type column chart for the allele frequencies of the insetions
g2 <- ggplot(Merge_sort4,
            aes(x = freq,
                y = gene_name, fill=Sample_type)) + 
  geom_col(position = "identity", alpha=0.4) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() + 
  labs(title = "Allele frequency of insertion per gene in the case and control samples", x = "Frequency", y = "Gene name") +
  scale_x_continuous(breaks  = freq_range_seq,
                     labels = abs(round(freq_range_seq,0))) 
# Add a second scale to the y-axis to display the chr:pos information
g2 + guides(y.sec = guide_axis_manual(
  labels = levels(Merge_sort4$ID)))
