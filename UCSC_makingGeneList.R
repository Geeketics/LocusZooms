## UCSC gene-track export - 27-02-17
UCSC_GRCh37_Genes <- read.delim("~/GitHub/LocusZooms/UCSC_GRCh37_Genes.txt", stringsAsFactors=FALSE)

## exclude non-useful columns that were exported
keep_cols <- c("hg19.ensGene.chrom", "hg19.ensGene.txStart", "hg19.ensGene.txEnd", "hg19.ensGene.cdsStart", "hg19.ensGene.cdsEnd", "hg19.ensGene.strand", "hg19.ensGtp.gene", "hg19.ensGtp.transcript", "hg19.ensGtp.protein", "hg19.kgXref.geneSymbol")
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[,keep_cols]

## exclude ulternate scaffold chromosomes
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[grep(UCSC_GRCh37_Genes$hg19.ensGene.chrom, pattern = "_", invert = TRUE),]

## exclude directly duplicated lines
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[!duplicated(UCSC_GRCh37_Genes),]

## add extra gene info to data.frame
UCSC_GRCh37_Genes$gene.length = abs(UCSC_GRCh37_Genes$hg19.ensGene.txEnd - UCSC_GRCh37_Genes$hg19.ensGene.txStart)
UCSC_GRCh37_Genes$cds.length = abs(UCSC_GRCh37_Genes$hg19.ensGene.cdsEnd - UCSC_GRCh37_Genes$hg19.ensGene.cdsStart)
UCSC_GRCh37_Genes$Coding[UCSC_GRCh37_Genes$cds.length == 0] = "Non-Coding"
UCSC_GRCh37_Genes$Coding[UCSC_GRCh37_Genes$cds.length > 0] = "Coding"

## determine which version of duplicate genes is longest
gene_names <- unique(UCSC_GRCh37_Genes$hg19.ensGtp.gene)

UCSC_GRCh37_Genes$Notes = NA
for(gene in gene_names){
  gene_set <- UCSC_GRCh37_Genes[UCSC_GRCh37_Genes$hg19.ensGtp.gene == gene,]
  longest <- max(gene_set$gene.length)
  UCSC_GRCh37_Genes$Notes[UCSC_GRCh37_Genes$hg19.ensGtp.gene == gene & UCSC_GRCh37_Genes$gene.length == longest] = "Longest Version"
  rm(gene_set, longest)
}

## keep longest version of each gene
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[!is.na(UCSC_GRCh37_Genes$Notes),]

## residual duplicates - just keep first instance
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[!duplicated(UCSC_GRCh37_Genes$hg19.ensGtp.gene),]

## rename columns for locuszoom.R
UCSC_GRCh37_Genes$Notes <- NULL
names(UCSC_GRCh37_Genes) <- c("Chrom", "Start", "End", "cdsStart", "cdsEnd", "Strand", "ensemblGeneID", "ensemblTranscriptID", "ensemblProteinID", "Gene", "GeneLength", "cdsLength", "Coding")

## change chromosome info to match PLINK convention (X = 23, Y = 24, MT = 26)
UCSC_GRCh37_Genes$Chrom <- gsub(UCSC_GRCh37_Genes$Chrom, pattern = "chr", replacement = "")
UCSC_GRCh37_Genes$Chrom[UCSC_GRCh37_Genes$Chrom == "X"] = 23
UCSC_GRCh37_Genes$Chrom[UCSC_GRCh37_Genes$Chrom == "Y"] = 24
UCSC_GRCh37_Genes$Chrom[UCSC_GRCh37_Genes$Chrom == "M"] = 26

## order based on start position & chromosome
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[order(UCSC_GRCh37_Genes$Start),]
UCSC_GRCh37_Genes <- UCSC_GRCh37_Genes[order(UCSC_GRCh37_Genes$Chrom),]

## save new file
write.table(UCSC_GRCh37_Genes, file = "~/GitHub/LocusZooms/UCSC_GRCh37_Genes_UniqueList.txt", quote = F, row.names = F, sep = "\t", na = "")

