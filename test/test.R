library(dplyr)
library(vroom)

Example.assoc.linear <- read.delim("Example.assoc.linear", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("Example.ld", stringsAsFactors = FALSE, header = TRUE)
UCSC_GRCh37_Genes_UniqueList.txt <- read.delim("UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = FALSE, header = TRUE)

# load the locuszoom function into R
source("functions/locus_zoom.R")

single = "rs2540781"
secondary = c("rs2540781", "rs8053279")

# create a LocusZoom-like plot
source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = NA, population = "EUR", sig.type = "P")

# Testing with one secondary variant
source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = single, population = "EUR", sig.type = "P")

# Testing with two secondary variant
source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs1008400", ld.file = Example.ld, offset = 200000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "Association of FTO with BMI in Europeans", nominal = 6, significant = 7.3, file.name = "Example.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P")

# Testing with a larger data set
Example.assoc.linear <- read.delim("eur_chr4.txt", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("rs145179124.ld", stringsAsFactors = FALSE, header = TRUE)

Example.assoc.linear$P = Example.assoc.linear$P + .Machine$double.xmin

source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs145179124", ld.file = Example.ld, offset = 500000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs145179124.jpg", secondary.snp = NA, population = "EUR", sig.type = "P")

# Try locating other loci as secondary SNP:
loci = read.delim('EUR_meta_full1_clean_rsid_0.01.clumped.clean', stringsAsFactors = F, header = T) %>% select(1:5)

Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs145179124')] = Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs145179124')] + .Machine$double.xmin
secondary = loci$SNP[which(loci$SNP != 'rs145179124')]

source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs145179124", ld.file = Example.ld, offset = 500000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs145179124.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P", secondary.label = T)

# Try with SLC2A9:
Example.assoc.linear <- read.delim("eur_chr4.txt", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("rs76242518.ld", stringsAsFactors = FALSE, header = TRUE)

Example.assoc.linear$P = Example.assoc.linear$P + .Machine$double.xmin

Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs76242518')] = Example.assoc.linear$P[which(Example.assoc.linear$SNP != 'rs76242518')] + .Machine$double.xmin
secondary = loci$SNP[which(loci$SNP != 'rs76242518')]

source("functions/locus_zoom.R")
locus.zoom(data = Example.assoc.linear, snp = "rs76242518", ld.file = Example.ld, offset = 500000, genes.data = UCSC_GRCh37_Genes_UniqueList.txt, noncoding = FALSE, plot.title = "EUR gout", nominal = 6, significant = 7.3, file.name = "rs76242518.jpg", secondary.snp = secondary, population = "EUR", sig.type = "P", secondary.label = T, nplots = 2)

