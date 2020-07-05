######################
# Locus Zoom, Make LD Data
# September 2016/February 2017
# Tanya Flynn
# Uni of Otago

### Important Running Notes:
## You need to have acces to the data that the associations were originally done on OR data to calculate ethnicity with
## Genes.Data expects a data.frame with headers Gene, Chrom, Start, & End - this will be used to annotate the graph.

## LD Colours
# 1.0 - 0.8 = #FF0000
# 0.8 - 0.6 = #FFA500
# 0.6 - 0.4 = #00FF00
# 0.4 - 0.2 = #87CEFA
# 0.2 - 0.0 = #000080
# NA = #7F7F7F
# top-hit = #7D26CD

# TODO: Add option to ignore lead variant, and use the specified SNP instead
locus.zoom = function(data = NULL, snp = NA, gene = NA, region = NA, ld.file = NULL, offset = 200000, genes.data = NULL, noncoding = FALSE, plot.title = NULL, nominal = 6, significant = 7.3, file.name = NULL, secondary.snp = NA, secondary.label = F, population = "EUR", sig.type = "P", nplots = 1)
{
	# Load Data and define constants:
	data$SNP = as.character(data$SNP)
	LD.colours = data.frame(LD = c(seq(from = 0, to = 1, by = 0.1), "NaN"), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2), "#7F7F7F"))

	# Stop if CHR:POS ID, but with no alias:
	if (!grepl('rs', snp) & is.null(alias)) {
		stop("You didn't provide an alias for your CHR:POS ID")
	}

	# Stop if there are duplicate SNPs:
	if (length(which(duplicated(data$SNP))) > 0) {
		stop('There are duplicates in your results file - Please remove them before running again')
	}

	# Get start and end regions for plotting and for pulling out data:
	if (all(is.na(region))) {
		# If SNP is given:
		if (!is.na(snp)) {
			snp.ind = which(data$SNP == snp)
			snp.chr = data$CHR[snp.ind]
			snp.pos = data$BP[snp.ind]
			region = c(snp.chr, snp.pos, snp.pos)
		} else if (!is.na(gene)) {
			# If Gene is given
			gene.ind = which(genes.data$Gene == gene)
			gene.chr = genes.data$Chrom[gene.ind]
			gene.start = genes.data$Start[gene.ind]
			gene.end = genes.data$End[gene.ind]
			region = c(gene.chr, gene.start, gene.end)
		} else {
			# If nothing was given
			stop("You must specify a SNP, Gene or Region to plot")
		}
	} else {
		offset = ifelse(is.na(offset), 0, offset)
	}

	# Now re-define region to work with:
	region[2] = region[2] - offset # start position
	region[3] = region[3] + offset # end position

	# Pull out the relevant information from the UCSC gene data.
	# Any gene that overlaps/intersect with the defined region is
	# included:
	genes.data = genes.data[genes.data$Chrom == region[1], ]
	genes.data = genes.data[genes.data$End > region[2], ]
	genes.data = genes.data[genes.data$Start < region[3], ]

	# Remove Non-Coding Gene Info:
	if(noncoding == FALSE){
		genes.data = genes.data[genes.data$Coding != "Non-Coding",]
	}

	# Likewise, pull out the relevant data from the results file:
	data = data[data$CHR == region[1], ]
	data = data[data$BP >= region[2], ]
	data = data[data$BP <= region[3], ]

	# Identify SNP with most significant association:
	lead.ind = which(data$P %in% min(data$P, na.rm = T))
	lead.snp = data$SNP[lead.ind]
	lead.chr = data$CHR[lead.ind]
	lead.pos = data$BP[lead.ind]
	lead.p = data$P[lead.ind]

	# If LD information is not supplied, calculate it from the 1000 genomes
	# data:
	if (is.null(ld.file)) {
		ld.file = get.ld(region, lead.snp, population)
	}

	# Check if results file is rsID-based or CHR:POS-based
	if (length(which(grepl('rs', data$SNP))) < 1) {
		ld.file$SNP_B = paste(ld.file$CHR_B, ld.file$BP_B, sep = ':')
	}

	# Add LD to Results
	new.ld.file = ld.file[,c("SNP_B", "R2")]

	round.up = function(x, decimals = 1){
		round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
	}

	data = merge(data, new.ld.file, by.x = "SNP", by.y = "SNP_B", all.x = TRUE)

	data$plot.ld = round.up(data$R2, decimals = 1)
	data$plot.ld[data$plot.ld > 1 & !is.na(data$plot.ld)] = 1

	data.plot = merge(data, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)

	# Make Plotting Variables
	y.max = max(-log10(lead.p), 8)
	x.min = region[2]
	x.max = region[3]

	# Make Plot

	# Define output plot size
	if (!is.integer(nplots) && nplots <= 0) {
		stop("You must specify number (whole integer) of loci to plot")
	}
	jpeg.height = (nplots * 80) + 40
	jpeg(width = 150, height = jpeg.height, units = "mm", res = 300, file = file.name)
	mat.row = (2 * nplots) + 1
	locus.par = c(4, 20)
	layout(matrix(c(1:mat.row), byrow = TRUE), heights = c(rep(locus.par, nplots), 8))

	# Plot N locus zooms
	# TODO: allow plotting of locus from different data sets
	plot.var = c(y.max, x.min, x.max, lead.snp, lead.pos, lead.p, nominal, significant)
	for (i in 1:nplots) {
		plot.locus(data.plot = data.plot, plot.title = plot.title, secondary.snp = secondary.snp, secondary.label = secondary.label, sig.type = sig.type, plot.var = plot.var)
	}

	# Plot Gene tracks
	par(mar = c(4, 4, 0.5, 4), mgp = c(2, 1, 0))
	plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", lead.chr), ylab="", xlim = c(x.min, x.max), ylim = c(0,2))

	# Stagger the genes
	y = rep(c(1.5, 0.5), times = length(genes.data$Gene))
	genes.data$Y = y[1:length(genes.data$Gene)]
	genes.top = genes.data[genes.data$Y == 1.5,]
	genes.bot = genes.data[genes.data$Y == 0.5,]

	# Plot the gene tracks:
	if (nrow(genes.top) > 0) {
		gene.position(genes.top)
	}
	if (nrow(genes.bot) > 0) {
		gene.position(genes.bot)
	}

	dev.off()
}

plot.locus <- function(data.plot = NULL, plot.title = NULL, nominal = 6, significant = 7.3, secondary.snp = NA, secondary.label = F, sig.type = "P", plot.var = NULL) {
	# Variables:
	# plot.var = c(y.max, x.min, x.max, lead.snp, lead.pos, lead.p, nominal, significant)
	y.max = as.numeric(plot.var[1])
	x.min = as.numeric(plot.var[2])
	x.max = as.numeric(plot.var[3])
	lead.snp = plot.var[4]
	lead.pos = as.numeric(plot.var[5])
	lead.p = as.numeric(plot.var[6])
	nominal = as.numeric(plot.var[7])
	significant = as.numeric(plot.var[8])

	# Plot SNP presence:
	par(mar = c(0, 4, 2, 4), mgp = c(2, 1, 0))
	plot(x = data.plot$BP, y = rep(1, times = nrow(data.plot)), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted\nSNPs", las = 2, xlim = c(x.min, x.max), cex.lab = 0.8)
	title(plot.title, line = 0)

	# Plot Manhattan/LocusZoom of region
	par(mar = c(0, 4, 0, 4), mgp = c(2, 1, 0))
	ylab = ifelse(sig.type == "P", expression(-log[10](italic(P))), expression(log[10](italic(BF))))
	plot(x = data.plot$BP, y = -log10(data.plot$P), ylim = c(0, y.max*1.1), pch = 20, col = as.character(data.plot$Colour), xlab = "", ylab = ylab, cex = 0.8, xaxt = "n", xlim = c(x.min, x.max))
	abline(h = nominal, col = "blue", lty = "dashed")
	abline(h = significant, col = "red", lty = "dashed")

	# Plot the lead SNP
	points(x = lead.pos, y = -log10(lead.p), pch = 18, cex = 2, col = "#7D26CD")
	text(x = lead.pos, y = -log10(lead.p), labels = lead.snp, pos = 3)

	# Plot label/text for the secondary SNP
	if(any(!is.na(secondary.snp))){
		check = which(data.plot$SNP %in% secondary.snp)
		if (length(check) != 0) {
			secondary.data = data.plot[check,]
			for (i in 1:nrow(secondary.data)) {
				plot.secondary.point(secondary.data, secondary.data$SNP[i], label = secondary.label)
			}
		} else {
			message('There was no secondary SNP in the region you wanted to plot:\n\tPlease check the location of your secondary SNP(s) - plotting without secondary SNP(s)')
			break
		}
	}

	# Add LD legend
	legend.colour = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080")
	legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2"), col = legend.colour, fill = legend.colour, border = legend.colour, pt.cex = 1.2, cex = 0.8, bg = "white", box.lwd = 0, title = expression("r"^2), inset = 0.01)
}

# Function to plot secondary SNP:
plot.secondary.point <- function (data, snp, label = F) {
	ind = which(data$SNP == snp)
	snp = data$SNP[ind]
	pos = data$BP[ind]
	logp = -log10(data$P[ind])
	points(x = pos, y = logp, pch = 1, cex = 1.1, col = "#FF0000")
	if (label) {
		text(x = pos, y = logp, labels = snp, cex = 0.7, pos = 3)
	}
}

# Function to plot the gene tracks and labels properly:
# TODO: double-check this works properly (e.g. no overlapping labels)
gene.position = function(data) {
	odd = 1
	for (i in 1:length(data$Gene)) {
		lines(x = c(data$Start[i], data$End[i]), y = c(data$Y[i], data$Y[i]), lwd = 3, col = "#000080")
		if(length(data$Gene) >= 10){
			length = abs(data$Start[i] - data$End[i])
			if(length > 8000) {
				if (odd%%2 == 0) {
					text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] - 0.1, labels = data$Gene[i], cex = 0.8, pos = 3, cex = 0.8)
				} else {
					# TODO: double-check the pos argument
					text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i], labels = data$Gene[i], cex = 0.8, pos = 1, cex = 0.8)
				}
			}
		} else {
			text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] - 0.1, labels = data$Gene[i], pos = 3, cex = 0.8)
		}
		odd = odd + 1
	}
}

# Function to get the LD information of specified population from the 1000
# Genomes data (March 2017 release).
#
# Note that the input SNP MUST be in rsID format, not CHR:POS-based.
#
# This function will leave/save the LD information in the working directory for
# future reference (e.g. if the user wanted to use the same LD information)
get.ld = function(region, snp, population) {
ld.snp = snp
# If the SNP is in CHR:POS ID, then find the rsID:
if (!grepl('rs', snp)) {
	command = "awk '{print $1\":\"$2, $3}' /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/1kgp_chrZZ_biallelic_snps.txt | grep -w SNP"
	command = gsub('ZZ', region[1], command)
	command = gsub('SNP', snp, command)
	rsid = system(command, intern = T)
	rsid = unlist(strsplit(rsid, ' '))[2]
	ld.snp = rsid
}

	# gsub the command and filename for chr, start/end positions and the
	# population:
	base.command = "bcftools view -r ZZ:Y1-Y2 -S /Volumes/archive/merrimanlab/riku/1kgp_sample_files/POP.sample /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/ALL.chrZZ.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o tmp.vcf.gz && plink2 --vcf tmp.vcf.gz --allow-no-sex --snps-only --r2 --inter-chr --ld-snp SNP --ld-window-r2 0 --out POP_region_ZZ_Y1_Y2 && rm tmp.vcf.gz POP_region_ZZ_Y1_Y2.nosex"
	base.command = gsub('ZZ', region[1], base.command)
	base.command = gsub('Y1', region[2], base.command)
	base.command = gsub('Y2', region[3], base.command)
	base.command = gsub('POP', population, base.command)
	base.command = gsub('SNP', ld.snp, base.command)

	# Make a system call to run the bcftools/plink command.
	# I'm only assigning it to a variable to suppress any possible form of
	# messages/outputs from the command, just in case
	messages = system(base.command, ignore.stdout = T, intern = T)

	# Import the LD data:
	ld.file = "POP_region_ZZ_Y1_Y2.ld"
	ld.file = gsub('ZZ', region[1], ld.file)
	ld.file = gsub('Y1', region[2], ld.file)
	ld.file = gsub('Y2', region[3], ld.file)
	ld.file = gsub('POP', population, ld.file)

	# Import the ld file and return the data
	ld = read.table(ld.file, stringsAsFactors = F, header = T)
	return(ld)
}
