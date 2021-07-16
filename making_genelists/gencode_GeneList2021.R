## gencode annotation download - 25-06-21
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz

library(tidyverse)

gencode.v38lift37.annotation <- read.delim("making_genelists/gencode.v38lift37.annotation.gtf")
  # original: gencode.v38lift37.annotation <- read.delim("~/Downloads/gencode.v38lift37.annotation.gtf", header=FALSE, comment.char="#", quote = "")
# column names acquired from: https://www.gencodegenes.org/pages/data_format.html
names(gencode.v38lift37.annotation) <- c("chromosome_name", "annotation_source", "feature_type", "genomic_start", "genomic_end", "score", "genomic_strand", "genomic_phase", "tosplit")


#### convert final column into multiple columns (per annotation) ####
# add index column to track original row of info
gencode.v38lift37.annotation$index <- 1:length(gencode.v38lift37.annotation[, 1])

# split list of annotations - new annotation per row (into label & value as a single column)
gene_IDs <- gencode.v38lift37.annotation %>%
  separate_rows(tosplit, sep = "; ")

# save first stage file
write.table(gene_IDs, "temp_gencode.annotation.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE, col.names = FALSE)
  # original: write.table(gene_IDs, "~/Downloads/gencode.v38lift37.annotation.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE, col.names = FALSE)

# split label and value into new columns with bash
system("sed 's/ /\t/g' temp_gencode.annotation.txt > temp_gencode.annotation.split.txt")
  # original: sed 's/ /\t/g' ~/Downloads/gencode.v38lift37.annotation.txt > ~/Downloads/gencode.v38lift37.annotation.split.txt

# load edited file
gene_IDs <- read.delim("temp_gencode.annotation.split.txt", header=FALSE)
  # original: gene_IDs <- read.delim("~/Downloads/gencode.v38lift37.annotation.split.txt", header=FALSE)

names(gene_IDs) <- c("chromosome_name", "annotation_source", "feature_type", "genomic_start", "genomic_end", "score", "genomic_strand", "genomic_phase", "label", "value", "index")

# convert to wide-format, using index column to get back to original number of rows
gene_IDs2 <- gene_IDs %>%
  pivot_wider(id_cols = index, names_from = label, values_from = value, values_fn = list)

# convert out of tibble format, remove lists of values
gene_IDs3 <- as.data.frame(gene_IDs2)
for(col in 1:length(gene_IDs3[1, ])){
  gene_IDs3[, col] <- sapply(gene_IDs3[, col], function(x) paste(x, collapse = "; "))
}
rm(col)

# add original info columns into new wide format file
gencode.wide <- merge(gencode.v38lift37.annotation[, -9], gene_IDs3, by = "index")

# save converted version of file
write.table(gencode.wide, "making_genelists/gencode.v38lift37.annotation.wide.gtf", quote=FALSE, sep = "\t", na = "", row.names = FALSE)
  # original: write.table(gencode.wide, "~/Downloads/gencode.v38lift37.annotation.wide.gtf", quote=FALSE, sep = "\t", na = "", row.names = FALSE)

# remove intermediate steps
rm(gene_IDs, gene_IDs2, gene_IDs3)
system("rm temp_gencode*")
# original: rm ~/Downloads/gencode.v38lift37.*.txt


#### create unique gene list for locuszooms etc ####
# exclude non-useful columns
keep_cols <- c("chromosome_name", "genomic_start", "genomic_end", "genomic_strand", "gene_id", "transcript_id", "protein_id", "gene_name", "gene_type", "transcript_type", "feature_type")

gencode.wide.working <- gencode.wide[, keep_cols]

rm("keep_cols")

# exclude directly duplicated lines
gencode.wide.working <- gencode.wide.working[!duplicated(gencode.wide.working), ]

# keep lines for genes, transcripts, and CDS only
gencode.wide.working <- gencode.wide.working[gencode.wide.working$feature_type %in% c("gene", "transcript", "CDS"), ]

# relabel gene_type categories (definitions from: https://www.gencodegenes.org/pages/biotypes.html)
gene_type_proteincoding <- c("protein_coding", "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene", "TR_D_gene", "TEC")
gene_type_psuedogene <- c("transcribed_unprocessed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "pseudogene", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "unitary_pseudogene", "translated_processed_pseudogene", "IG_C_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "IG_J_pseudogene", "polymorphic_pseudogene", "rRNA_pseudogene")
gene_type_longnoncodingRNA <- c("lncRNA",  "lincRNA", "sense_overlapping", "sense_intronic", "antisense", "processed_transcript")
gene_type_noncodingRNA <- c("snRNA", "miRNA", "misc_RNA", "rRNA","scRNA", "Mt_rRNA", "Mt_tRNA", "snoRNA", "vault_RNA")

gencode.wide.working$gene_type[gencode.wide.working$gene_type %in% gene_type_proteincoding] <- "proteincoding"
gencode.wide.working$gene_type[gencode.wide.working$gene_type %in% gene_type_psuedogene] <- "psuedogene"
gencode.wide.working$gene_type[gencode.wide.working$gene_type %in% gene_type_longnoncodingRNA] <- "lncRNA"
gencode.wide.working$gene_type[gencode.wide.working$gene_type %in% gene_type_noncodingRNA] <- "ncRNA"

rm("gene_type_longnoncodingRNA", "gene_type_noncodingRNA", "gene_type_proteincoding", "gene_type_psuedogene")

# convert each gene_id into a single line
gencode.wide.unique <- gencode.wide.working[1, ]
gencode.wide.unique[, c("genomic_length", "cds_start", "cds_end", "cds_length", "coding")] <- NA
gencode.wide.unique[, c("gene_type", "transcript_type", "feature_type" )] <- NULL
gencode.wide.unique <- gencode.wide.unique[0, ]

row <- 1
for (gene in unique(gencode.wide.working$gene_id)){
  info <- gencode.wide.working[gencode.wide.working$gene_id == gene, ]
  info[info == ""] <- NA
  
  # chromosome
  gencode.wide.unique[row, "chromosome_name"] <- paste(unique(info$chromosome_name[!is.na(info$chromosome_name)]), collapse = "; ")
  # start
  gencode.wide.unique[row, "genomic_start"] <- min(info$genomic_start, na.rm = TRUE)
  # end
  gencode.wide.unique[row, "genomic_end"] <- max(info$genomic_end, na.rm = TRUE)
  # CDSstart
  if(length(info[info$feature_type == "CDS", 1]) > 0) { gencode.wide.unique[row, "cds_start"] <- min(info$genomic_start[info$feature_type == "CDS"], na.rm = TRUE) }
  # CDSend
  if(length(info[info$feature_type == "CDS", 1]) > 0) { gencode.wide.unique[row, "cds_end"] <- max(info$genomic_end[info$feature_type == "CDS"], na.rm = TRUE) }
  # strand
  gencode.wide.unique[row, "genomic_strand"] <- paste(unique(info$genomic_strand[!is.na(info$genomic_strand)]), collapse = "; ")
  # gene_id
  gencode.wide.unique[row, "gene_id"] <- paste(unique(info$gene_id[!is.na(info$gene_id)]), collapse = "; ")
  # transcript_id
  gencode.wide.unique[row, "transcript_id"] <- paste(unique(info$transcript_id[!is.na(info$transcript_id)]), collapse = "; ")
  # protein_id
  gencode.wide.unique[row, "protein_id"] <- paste(unique(info$protein_id[!is.na(info$protein_id)]), collapse = "; ")
  # gene_name
  gencode.wide.unique[row, "gene_name"] <- paste(unique(info$gene_name[!is.na(info$gene_name)]), collapse = "; ")
  # max length
  gencode.wide.unique[row, "genomic_length"] <- abs(gencode.wide.unique[row, "genomic_end"] - gencode.wide.unique[row, "genomic_start"])
  # CDS length
  gencode.wide.unique[row, "cds_length"] <- abs(gencode.wide.unique[row, "cds_end"] - gencode.wide.unique[row, "cds_start"])
  # coding
  gencode.wide.unique[row, "coding"] <- paste(unique(info$gene_type[!is.na(info$gene_type)]), collapse = "; ")
  row <- row + 1
  rm(info)
}
rm(row, gene)

# match column header to UCSC sourced list (from LocusZoom code)
ucsc.unique <- read.delim("UCSC_GRCh37_Genes_UniqueList2017.txt", header = TRUE, stringsAsFactors = FALSE)
  # original: ucsc.unique <- read.delim("~/Documents/GitHubs/LocusZooms/UCSC_GRCh37_Genes_UniqueList.txt", header = TRUE, stringsAsFactors = FALSE)

names(ucsc.unique)

gencode.wide.unique <- gencode.wide.unique[, c(1:3, 10:11, 4:9, 12:13)]
names(gencode.wide.unique) <- names(ucsc.unique)

## order based on start position & chromosome
gencode.wide.unique <- gencode.wide.unique[order(factor(gencode.wide.unique$Chrom, levels = paste0("chr", c(1:22, "X", "Y", "M"))), gencode.wide.unique$Start), ]

# save file of unique genes
write.table(gencode.wide.unique, "Gencode_GRCh37_Genes_UniqueList2021.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE)
  # original: write.table(gencode.wide.unique, "~/Documents/GitHubs/LocusZooms/Gencode_GRCh37_Genes_UniqueList.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE)

rm("gencode.wide", "gencode.wide.working")


