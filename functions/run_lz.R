# Rscript to generate locus zooms for all the loci listed in PLINK clump output

library(vroom)
library(dplyr)
library(furrr)

plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3)

# Load LZ functions
source('src/lz_scripts/locus_zoom.R')

gene_list = read.delim("src/lz_scripts/UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = F, header = T)

args = commandArgs(trailingOnly = T)

loci = read.plink.loci(args[1])

ancestry = gsub('results/locus_zooms/loci/', '', args[1])
ancestry = gsub('/.*', '', ancestry)

# Load in data:
if (length(args) > 2) {
	data = list()
	for (i in 2:length(args)) {
		list[[i - 1]] = vroom(args[i], col_types = cols(P = col_character()))
	}
	names(data) = args[2:length(args)]
} else {
	data = vroom(args[2], col_types = cols(P = col_character()))
}

# TODO: make sure urate data is clean before loading
# urate = read.delim('eur_chr4_urate.txt', stringsAsFactors = F, header = T, sep = ' ')
# colnames(urate) = c( "CHR", "BP", "SNP", "Allele1", "Allele2", "MAF", "Effect", "StdErr", "P", "N")

# Load all the LD info:
ld_path = paste('results/locus_zooms/ld/', ancestry, sep = '')
ld_path = paste(ld_path, '/', sep = '')

ld_files = paste(ld_path, loci$SNP, sep = '')
ld_files = paste(ld_files, '.ld', sep = '')

ld_list = lapply(ld_files, function(x) read.table(x, stringsAsFactors = F, header = T))

names(ld_list) = loci$SNP

secondary = loci$SNP

out_dir = paste('results/locus_zooms/lz_plots/', ancestry, sep = '')

plot_query <- function(data, ld_list, gene_list, query_snp, secondary, ancestry, out_dir) {
	if (is.data.frame(data)) {
		plot_title = paste(c(ancestry, "GOUT", query_snp), collapse = ' ')
		nplots = F
	} else {
		plot_title = paste(ancestry, names(data), sep = ' ')
		plot_title = paste(plot_title, query_snp, sep = ' ')
		nplots = T
	}
	jpeg_name = paste(query_snp, '.jpg', sep = '')
	jpeg_name = paste(out_dir, jpeg_name, sep = '/')
	locus.zoom(data = data, snp = query_snp, ld.file = ld_list[[query_snp]], offset = 500000, genes.data = gene_list, noncoding = FALSE, plot.title = plot_title, nominal = 6, significant = 7.3, file.name = jpeg_name, secondary.snp = secondary, population = ancestry, sig.type = "P", secondary.label = T, nplots = nplots, ignore.lead = T, rsid.check = F)
}

# Plot locus zooms if the number of loci to plot is > 50:
# (With less loci, I think it will run faster with for-loop, as furrr takes
# a while to split/initialise the jobs into multiple processes)
if (length(loci$SNP) > 50) {
	future_map(loci$SNP, function(x) plot_query(data, ld_list, gene_list, x, secondary, ancestry, out_dir))
} else {
	for (i in 1:length(loci$SNP)) {
		query_snp = loci$SNP[i]
		plot_query(data, ld_list, gene_list, query_snp, secondary, ancestry, out_dir)
	}
}

