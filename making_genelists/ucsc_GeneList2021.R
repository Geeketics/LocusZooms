#### create unique gene list for locuszooms etc ####
# exported from UCSC gene tracks - 5 Jul 2021
big_ucsc <- read.delim("making_genelists/ucsc_ensembl_export_50721.txt")
  # original: big_ucsc <- read.delim("~/Downloads/ucsc_ensembl_export (2)")

# exclude non-useful columns
keep_cols <- c("hg19.ensGene.chrom", "hg19.ensGene.txStart", "hg19.ensGene.txEnd", "hg19.ensGene.cdsStart", "hg19.ensGene.cdsEnd", "hg19.ensGene.strand", "hg19.ensGtp.gene", "hg19.ensGtp.transcript", "hg19.ensGtp.protein", "hg19.ensemblToGeneName.value", "hg19.ensemblSource.source")

big_ucsc.working <- big_ucsc[, keep_cols]

rm("keep_cols")

# exclude directly duplicated lines
big_ucsc.working <- big_ucsc.working[!duplicated(big_ucsc.working), ]

# relabel ensemble source categories (definitions from: https://www.gencodegenes.org/pages/biotypes.html)
gene_type_proteincoding <- c("protein_coding", "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_V_gene", "TR_C_gene", "TR_V_gene", "TR_J_gene", "TR_D_gene", "non_stop_decay", "nonsense_mediated_decay")
gene_type_psuedogene <- c("transcribed_unprocessed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "pseudogene", "transcribed_processed_pseudogene", "unitary_pseudogene", "translated_processed_pseudogene", "IG_C_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "IG_J_pseudogene", "polymorphic_pseudogene")
gene_type_longnoncodingRNA <- c("lincRNA", "sense_overlapping", "sense_intronic", "antisense",  "3prime_overlapping_ncrna", "retained_intron", "processed_transcript")
gene_type_noncodingRNA <- c("snRNA", "miRNA", "misc_RNA", "rRNA", "Mt_rRNA", "Mt_tRNA", "snoRNA")

big_ucsc.working$hg19.ensemblSource.source[big_ucsc.working$hg19.ensemblSource.source %in% gene_type_proteincoding] <- "proteincoding"
big_ucsc.working$hg19.ensemblSource.source[big_ucsc.working$hg19.ensemblSource.source %in% gene_type_psuedogene] <- "psuedogene"
big_ucsc.working$hg19.ensemblSource.source[big_ucsc.working$hg19.ensemblSource.source %in% gene_type_longnoncodingRNA] <- "lncRNA"
big_ucsc.working$hg19.ensemblSource.source[big_ucsc.working$hg19.ensemblSource.source %in% gene_type_noncodingRNA] <- "ncRNA"

rm("gene_type_longnoncodingRNA", "gene_type_noncodingRNA", "gene_type_proteincoding", "gene_type_psuedogene")

# convert each gene_id into a single line
big_ucsc.unique <- big_ucsc.working[1, ]
big_ucsc.unique[, c("genomic_length", "cds_length")] <- NA
big_ucsc.unique <- big_ucsc.unique[0, ]

row <- 1
for (gene in unique(big_ucsc.working$hg19.ensGtp.gene)){
  info <- big_ucsc.working[big_ucsc.working$hg19.ensGtp.gene == gene, ]
  info[info == ""] <- NA
  
  # chromosome
  big_ucsc.unique[row, "hg19.ensGene.chrom"] <- paste(sort(unique(info$hg19.ensGene.chrom[!is.na(info$hg19.ensGene.chrom)])), collapse = "; ")
  # start
  big_ucsc.unique[row, "hg19.ensGene.txStart"] <- min(info$hg19.ensGene.txStart, na.rm = TRUE)
  # end
  big_ucsc.unique[row, "hg19.ensGene.txEnd"] <- max(info$hg19.ensGene.txEnd, na.rm = TRUE)
  # CDSstart
  big_ucsc.unique[row, "hg19.ensGene.cdsStart"] <- min(info$hg19.ensGene.cdsStart, na.rm = TRUE)
  # CDSend
  big_ucsc.unique[row, "hg19.ensGene.cdsEnd"] <- max(info$hg19.ensGene.cdsEnd, na.rm = TRUE)
  # strand
  big_ucsc.unique[row, "hg19.ensGene.strand"] <- paste(sort(unique(info$hg19.ensGene.strand[!is.na(info$hg19.ensGene.strand)])), collapse = "; ")
  # gene_id
  big_ucsc.unique[row, "hg19.ensGtp.gene"] <- paste(sort(unique(info$hg19.ensGtp.gene[!is.na(info$hg19.ensGtp.gene)])), collapse = "; ")
  # transcript_id
  big_ucsc.unique[row, "hg19.ensGtp.transcript"] <- paste(sort(unique(info$hg19.ensGtp.transcript[!is.na(info$hg19.ensGtp.transcript)])), collapse = "; ")
  # protein_id
  big_ucsc.unique[row, "hg19.ensGtp.protein"] <- paste(sort(unique(info$hg19.ensGtp.protein[!is.na(info$hg19.ensGtp.protein)])), collapse = "; ")
  # gene_name
  big_ucsc.unique[row, "hg19.ensemblToGeneName.value"] <- paste(sort(unique(info$hg19.ensemblToGeneName.value[!is.na(info$hg19.ensemblToGeneName.value)])), collapse = "; ")
  # max length
  big_ucsc.unique[row, "genomic_length"] <- abs(big_ucsc.unique[row, "hg19.ensGene.txEnd"] - big_ucsc.unique[row, "hg19.ensGene.txStart"])
  # CDS length
  big_ucsc.unique[row, "cds_length"] <- abs(big_ucsc.unique[row, "hg19.ensGene.cdsEnd"] - big_ucsc.unique[row, "hg19.ensGene.cdsStart"])
  # coding
  big_ucsc.unique[row, "hg19.ensemblSource.source"] <- paste(sort(unique(info$hg19.ensemblSource.source[!is.na(info$hg19.ensemblSource.source)])), collapse = "; ")
  row <- row + 1
  rm(info)
}
rm(row, gene)

# remove weird chromosome scaffolds
big_ucsc.unique <- big_ucsc.unique[!grepl(big_ucsc.unique$hg19.ensGene.chrom, pattern = "_"), ]

# match column header to original UCSC sourced list (from LocusZoom code)
ucsc.unique <- read.delim("UCSC_GRCh37_Genes_UniqueList2017.txt", header = TRUE, stringsAsFactors = FALSE)
  # original: ucsc.unique <- read.delim("~/Documents/GitHubs/LocusZooms/UCSC_GRCh37_Genes_UniqueList.txt", header = TRUE, stringsAsFactors = FALSE)

names(ucsc.unique)

big_ucsc.unique <- big_ucsc.unique[, c(1:10, 12:13, 11)]
names(big_ucsc.unique) <- names(ucsc.unique)

# clarify gene coding classification where more than one was found (need one for LZ gene filtering code)
# base on hierarchy of classification 'proteincoding' > 'psuedogene' > 'lncRNA' > 'ncRNA'
big_ucsc.unique$Coding[grepl(big_ucsc.unique$Coding, pattern = "proteincoding")] <- "proteincoding"
big_ucsc.unique$Coding[grepl(big_ucsc.unique$Coding, pattern = "psuedogene")] <- "psuedogene"
big_ucsc.unique$Coding[grepl(big_ucsc.unique$Coding, pattern = "lncRNA")] <- "lncRNA"

## order based on start position & chromosome
big_ucsc.unique <- big_ucsc.unique[order(factor(big_ucsc.unique$Chrom, levels = paste0("chr", c(1:22, "X", "Y", "M"))), big_ucsc.unique$Start), ]

# save file of unique genes
write.table(big_ucsc.unique, "UCSC_GRCh37_Genes_UniqueList2021.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE)
  # original: write.table(big_ucsc.unique, "~/Documents/GitHubs/LocusZooms/UCSC_GRCh37_Genes_UniqueList2021.txt", quote = FALSE, sep = "\t", na = "", row.names = FALSE)


rm("big_ucsc", "big_ucsc.working", "ucsc.unique")


