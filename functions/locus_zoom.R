######################
# Locus Zoom, Make LD Data
# July 2020
# Tanya Major & Riku Takei
# Uni of Otago

### Important Running Notes:
## data: expects a data.frame (or a list of data.frames) with at least columns containing the chromosome, positions, rsIDs, and p-values of your results - these must be labelled CHR, BP, SNP, and P
## snp: specifies a SNP to centre the plot around
## gene: specifies a gene to centre the plot around
## region: specifies the chromosome start and end you wish to plot
## ld.file: expects a data.frame of the LD between your lead SNP and all other SNPs (requires the columns SNP_B and R2) - if left blank the script can calculate this for you so long as you have access to the biochem servers
## offset_bp: specifies how far either side of your gene/snp/region of interest to plot (in base pairs)
## genes.data: specifies a data.frame of genes within the plot region - expects the headers Gene, Chrom, Start, & End
## non-coding: specifies whether to annotate non-coding genes under the plot
## plot.title: specifies what to label the plot as
## nominal: specifies where to draw a nominal significance line
## significant: specifies where to draw a significance line
## file.name: specifies what name to save the plot under
## secondary.snp: specifies a list of SNPs to label on the graph as well as labelling the top SNP / chosen SNP
## secondary.label: specifies whether to label the secondary SNPs on the plot
## genes.pvalue: specifies a data.frame of p-values associated with each gene (e.g. MAGMA results) - expects the headers Gene, P
## colour.genes: specifies whether to colour genes based on a p-value provided in gene.pvalue
## population: specifies the 1000 genomes population to use for LD (if script is calculating for you)
## sig.type: specifies whether the y-axis should be -log10(P) or -log10(BF) - these are the only two options
## nplots: specifies how many plots will be saved into a single jpeg (e.g. plot two GWAS results one above another, nplots = TRUE)
## ignore.lead: specifies whether to ignore the SNP with the smallest P and use the SNP specified by 'snp' to centre the plot
## rsid.check: specifies whether to check if the SNPs are labelled with rsIDs - should only matter if script is calculating LD for you



## LD Colours
# 1.0 - 0.8 = #FF0000
# 0.8 - 0.6 = #FFA500
# 0.6 - 0.4 = #00FF00
# 0.4 - 0.2 = #87CEFA
# 0.2 - 0.0 = #000080
# NA = #7F7F7F
# top-hit = #7D26CD

## Gene Colours
# <1e-15 = #FF0000
# ≥1e-15, <1e-10 = #FFA500
# ≥1e-10, <1e-5 = #00FF00
# ≥1e-5, <0.05 = #87CEFA
# ≥0.05 = #000080
# NA = #7F7F7F

# Function to make LocusZoom like plots
locus.zoom <- function(data = NULL, snp = NA, gene = NA, region = NA, ld.file = NULL, offset_bp = 200000, genes.data = NULL, noncoding = FALSE, plot.title = NULL, plot.type = "jpg", nominal = 6, significant = 7.3, file.name = NULL, secondary.snp = NA, secondary.label = FALSE, genes.pvalue = NULL, colour.genes = FALSE, population = "EUR", sig.type = "P", nplots = FALSE, ignore.lead = FALSE, rsid.check = TRUE) {
  
  # Define constants:
  LD.colours <- data.frame(LD = as.character(seq(from = 0, to = 1, by = 0.1)), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2)), stringsAsFactors = FALSE)
  GENE.colours <- data.frame(Threshold = c(">2.6e-6", ">1e-10", ">1e-15", ">1e-20", "<1e-20"), Colour = c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), stringsAsFactors = FALSE)
  
  # Load Data
  # If plotting multiple summary stats, take the first summary stats as lead/reference data:
  if (is.data.frame(data)) {
    lead.data = data
  } else {
    lead.data = data[[1]]
  }
  lead.data$SNP = as.character(lead.data$SNP)
  
  # Check the SNPs are in rsID format and no duplicates:
  if (rsid.check) {
    check.rsid(lead.data$SNP)
  }
  
  # Get start and end regions for plotting and for pulling out data:
  if (all(is.na(region))) {
    region = get.region(lead.data, snp, genes.data, gene)
  } else {
    offset_bp = ifelse(is.na(offset_bp), 0, offset_bp)
  }
  
  # Now re-define region to work with:
  region[2] = region[2] - offset_bp # start position
  region[3] = region[3] + offset_bp # end position
  
  # Pull out the relevant information from the UCSC gene data.
  # Any gene that overlaps/intersect with the defined region is included:
  genes.data = genes.data[genes.data$Chrom == region[1], ]
  genes.data = genes.data[genes.data$End > region[2], ]
  genes.data = genes.data[genes.data$Start < region[3], ]
  
  # Remove Non-Coding Gene Info:
  if(!noncoding) {
    genes.data = genes.data[genes.data$Coding != "Non-Coding", ]
  }
  
  # Pull out the relevant information from the gene p-values data
  if(colour.genes) {
    genes.pvalue = genes.pvalue[genes.pvalue$Gene %in% genes.data$Gene, ]
  }
    
  # Pull out the relevant data from the result file(s), and logBF/log p-value:
  if (is.data.frame(data)) {
    data = subset.data(data, region)
    if (sig.type == "P") {
      data$logP = as.numeric(unlist(lapply(data$P, elog10)))
    } else {
      data$logP = data$logBF
    }
    lead.data = data
  } else {
    data = lapply(data, function(x) subset.data(x, region))
    lead.data = data[[1]]
  }
  
  # Get info on lead variant:
  if (ignore.lead) {
    if (is.na(snp)) {
      stop("You must provide a SNP with the ignore.lead = TRUE option")
    }
    lead.ind = which(lead.data$SNP == snp)
  } else {
    lead.ind = which(lead.data$logP %in% max(lead.data$logP, na.rm = TRUE))[1]
  }
  
  lead.snp = lead.data$SNP[lead.ind]
  lead.chr = lead.data$CHR[lead.ind]
  lead.pos = lead.data$BP[lead.ind]
  lead.logp = lead.data$logP[lead.ind]
  
  # If LD information is not supplied, calculate it from the 1000 genomes data:
  if (is.null(ld.file)) {
    ld.file = get.ld(region, lead.snp, population)
  }

  # Add LD to Results
  new.ld.file = ld.file[, c("SNP_B", "R2")]
  
  if (is.data.frame(data)) {
    data.plot = merge.plot.dat(data, new.ld.file, LD.colours)
  } else {
    data.plot = lapply(data, function(x) merge.plot.dat(x, new.ld.file, LD.colours))
  }
  
  # Make Plot
  
  # Define output plot size
  npanel = ifelse(nplots, length(data), 1)
  plot.height = (npanel * 80) + 50
  if(plot.type == "jpg"){
    jpeg(width = 160, height = plot.height, units = "mm", res = 300, filename = file.name)
  } else{
    svg(width = (160 / 25.4), height = (plot.height / 25.4), filename = file.name)
  }

  mat.row = (2 * npanel) + 1
  locus.par = c(4, 20)
  layout(matrix(c(1:mat.row), byrow = TRUE), heights = c(rep(locus.par, npanel), 10))
  
  # Make Plotting Variables
  x.min = region[2]
  x.max = region[3]
  
  # Plot N locus zooms
  if (npanel > 1) {
    for (i in 1:npanel) {
      # Set y.max:
      tmp.dat = data.plot[[i]]
      y.max = max(tmp.dat$logP, 8)
      plot.var = c(y.max, x.min, x.max, lead.snp, nominal, significant)
      plot.locus(data.plot = tmp.dat, plot.title = names(data.plot)[i], secondary.snp = secondary.snp, secondary.label = secondary.label, sig.type = sig.type, plot.var = plot.var)
      rm(tmp.dat)
    }
  } else {
    y.max = max(data.plot$logP, 8)
    plot.var = c(y.max, x.min, x.max, lead.snp, nominal, significant)
    plot.locus(data.plot = data.plot, plot.title = plot.title, secondary.snp = secondary.snp, secondary.label = secondary.label, sig.type = sig.type, plot.var = plot.var)
  }
  
  # Plot Gene tracks
  par(mar = c(4, 4, 0.5, 8), mgp = c(2, 1, 0), xpd = FALSE)
  plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", lead.chr), ylab="", xlim = c(x.min, x.max), ylim = c(0, 3))

  if (nrow(genes.data) != 0) {
                       
  # add colour column to genes.data
  if(colour.genes) {
    genes.data = merge.gene.colour(genes.data, genes.pvalue, GENE.colours)
  } else {
    genes.data$Colour = "#7F7F7F"
  }
  
  # Stagger the genes
    y = rep(c(2.5, 1.5, 0.5), times = length(genes.data[ ,"Gene"]))
    genes.data$Y = y[1:length(genes.data$Gene)]
    genes.top <- genes.data[genes.data$Y == 2.5, ]
    genes.mid <- genes.data[genes.data$Y == 1.5, ]
    genes.bot <- genes.data[genes.data$Y == 0.5, ]
    
    # Plot the gene tracks:
    if (nrow(genes.top) > 0) {
      gene.position(genes.top)
    }
    if (nrow(genes.mid) > 0) {
      gene.position(genes.mid)
    }
    if (nrow(genes.bot) > 0) {
      gene.position(genes.bot)
    }
  }

  # add gene colour legend
  if(colour.genes) {
    legend.colour = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080", "#7F7F7F")
    par(xpd = TRUE)
    legend(x = "right", legend = c(expression("<1x10"^-20), expression(paste("<1x10"^-15, "; ≥1x10"^-20)), expression(paste("<1x10"^-10, "; ≥1x10"^-15)), expression(paste("<2.6x10"^-6, "; ≥1x10"^-10)), expression("≥2.6x10"^-6), "Unknown"), col = legend.colour, fill = legend.colour, border = legend.colour, pt.cex = 1.2, cex = 0.8, bg = "white", box.lwd = 0, title = "p-value", inset = -0.22)
  }

  dev.off()
}

# Function to create LocusZoom style plot (without gene track):
plot.locus <- function(data.plot = NULL, plot.title = NULL, nominal = 6, significant = 7.3, secondary.snp = NA, secondary.label = FALSE, sig.type = "P", plot.var = NULL, ignore.lead = FALSE) {
  # Variables:
  y.max = as.numeric(plot.var[1])
  x.min = as.numeric(plot.var[2])
  x.max = as.numeric(plot.var[3])
  lead.snp = plot.var[4]
  nominal = as.numeric(plot.var[5])
  significant = as.numeric(plot.var[6])
  
  # Plot SNP presence:
  par(mar = c(0, 4, 2, 8), mgp = c(2, 1, 0), xpd = FALSE)
  plot(x = data.plot$BP, y = rep(1, times = nrow(data.plot)), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted\nSNPs", las = 2, xlim = c(x.min, x.max), cex.lab = 0.8)
  title(plot.title, line = 0)
  
  # Plot Manhattan/LocusZoom of region
  par(mar = c(0, 4, 0, 8), mgp = c(2, 1, 0), xpd = FALSE)
  ylab = ifelse(sig.type == "P", expression(-log[10](italic(P))), expression(log[10](BF)))
  plot(x = data.plot$BP, y = data.plot$logP, ylim = c(0, y.max*1.1), pch = 20, col = as.character(data.plot$Colour), xlab = "", ylab = ylab, cex = 0.8, xaxt = "n", xlim = c(x.min, x.max))
  abline(h = nominal, col = "blue", lty = "dashed")
  abline(h = significant, col = "red", lty = "dashed")
  
  # Plot the lead SNP and remove it from the list of secondary SNPs, if present:
  if (lead.snp %in% data.plot$SNP) {
    ind = which(data.plot$SNP == lead.snp)
    lead.pos = data.plot$BP[ind]
    lead.logp = data.plot$logP[ind]
    points(x = lead.pos, y = lead.logp, pch = 18, cex = 1.5, col = "#7D26CD")
    text(x = lead.pos, y = lead.logp, labels = lead.snp, pos = 3)
    secondary.snp = secondary.snp[which(secondary.snp != lead.snp)]
  }
  
  # Plot label/text for the secondary SNP
  if(any(!is.na(secondary.snp))){
    check = which(data.plot$SNP %in% secondary.snp)
    if (length(check) != 0) {
      secondary.data = data.plot[check, ]
      for (i in 1:nrow(secondary.data)) {
        plot.secondary.point(data = secondary.data, snp = secondary.data$SNP[i], lead.snp = lead.snp, plot.var = plot.var, nominal = nominal, label = secondary.label)
      }
    }
  }
  
  # Add LD legend
  legend.colour = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080", "#7F7F7F")
  par(xpd = TRUE)
  legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2", "Unknown"), col = legend.colour, fill = legend.colour, border = legend.colour, pt.cex = 1.2, cex = 0.8, bg = "white", box.lwd = 0, title = expression("r"^2), inset = c(-0.14, 0.01))
}

# Function to check if the input variants have rsIDs
check.rsid <- function(snp = NULL) {
  # Stop if CHR:POS ID:
  if (all(!grepl('rs', snp))) {
    stop("Your SNP column does not have rsIDs")
  }
  # Stop if there are duplicate SNPs:
  if (length(which(duplicated(snp))) > 0) {
    stop("There are duplicate rsIDs in your results file - Please remove them before running again")
  }
}

# Function to make regions out of variant or gene information
get.region <- function(snp.dat, snp, gene.dat, gene) {
  # If SNP is given:
  if (!is.na(snp)) {
    snp.ind = which(snp.dat$SNP == snp)
    snp.chr = snp.dat$CHR[snp.ind]
    snp.pos = snp.dat$BP[snp.ind]
    region = c(snp.chr, snp.pos, snp.pos)
  } else if (!is.na(gene)) {
  # If Gene is given
    gene.ind = which(gene.dat$Gene == gene)
    gene.chr = gene.dat$Chrom[gene.ind]
    gene.start = gene.dat$Start[gene.ind]
    gene.end = gene.dat$End[gene.ind]
    region = c(gene.chr, gene.start, gene.end)
  } else {
  # If nothing was given
    stop("You must specify a SNP, Gene, or Region to plot")
  }
  return(region)
}

# Function to subset relevant region from a dataset
subset.data <- function(data, region) {
  res = data
  res = res[res$CHR == region[1], ]
  res = res[res$BP >= region[2], ]
  res = res[res$BP <= region[3], ]
  return(res)
}

round.up <- function(x, decimals = 1){
  round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
}

# Function to merge the LD and LD colours with the relevant region of the summary stats:
merge.plot.dat <- function(data, ld.file, LD.colours) {
  # add LD to data.frame with p-values
  res = merge(data, ld.file, by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
  # convert LD to categories
  res$plot.ld = round.up(res$R2, decimals = 1)
  res$plot.ld[res$plot.ld > 1 & !is.na(res$plot.ld)] = 1
  # add plotting colours based on LD categories
  res = merge(res, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)
  # make all variants without an LD value grey
  res$Colour[is.na(res$Colour)] = "#7F7F7F"
  # sort file based on position, then on LD (so LD > 0.2 not hidden by other points)
  res = res[order(res$BP), ]
  plot.last <- res[res$plot.ld > 0.2 & !is.na(res$plot.ld), ]
  res <- res[res$plot.ld <= 0.2 | is.na(res$plot.ld), ]
  res <- rbind(res, plot.last)
  rm(plot.last)
  return(res)
}

# Function to plot secondary SNP:
plot.secondary.point <- function(data, snp, lead.snp, plot.var, nominal, label = FALSE) {
  # save out variables
  ind = which(data$SNP == snp)
  snp = data$SNP[ind]
  pos = data$BP[ind]
  logp = data$logP[ind]
  lead.ind = which(data$SNP == lead.snp)
  lead.pos = data$BP[ind]

  # plot red line around secondary SNP
  points(x = pos, y = logp, pch = 1, cex = 1.1, col = "#FF0000")

  if (label) {
  # set up labelling offsets (x-axis)
    x.min = as.numeric(plot.var[2])
    x.max = as.numeric(plot.var[3])
    x.offset = abs(x.max - x.min) / 150 * 15
    
    if(pos < lead.pos) {
      label.x.offset = pos - x.offset
      line.x.offset = pos - (x.offset / 3)
      side = 1
    } else {
      label.x.offset = pos + x.offset
      line.x.offset = pos + (x.offset / 3)
      side = 0
    }
    
    # set up labelling offsets (y-axis) - considers SNPs around it  
    surrounding.data = data[data$BP > (pos - (abs(x.max - x.min) / 6)) & data$BP < (pos + (abs(x.max - x.min) / 6)), ]
    label.y.offset = max(surrounding.data$logP) * 1.03
    if(abs(nominal - label.y.offset) <= 1) {
      label.y.offset <- max(label.y.offset, nominal) * 1.03
    }
    
    # add lines and text to label SNP
    text(x = label.x.offset, y = (label.y.offset * 1.01), labels = snp, cex = 0.7, adj = c(1, side))
    segments(x0 = pos, x1 = line.x.offset, y0 = logp, y1 = label.y.offset)
  }
}


# Function to merge the gene and gene colours with the relevant region of the summary stats:
merge.gene.colour <- function(data, pvalues, GENE.colours) {
  # add gene pvalues to data.frame with gene positions
  res = merge(data, pvalues, by = "Gene", all.x = TRUE)
  # convert pvalues to categories
  res$gene.col = cut(res$P, breaks = c(0, 1e-20, 1e-15, 1e-10, 2.6e-6, 1), labels = c("<1e-20", ">1e-20", ">1e-15", ">1e-10", ">2.6e-6"), include.lowest = TRUE, right = TRUE)
  # add plotting colours based on pvalue categories
  res = merge(res, GENE.colours, by.x = "gene.col", by.y = "Threshold", all.x = TRUE)
  # make all genes without a pvalue grey
  res$Colour[is.na(res$Colour)] = "#7F7F7F"
  # sort by gene start position
  res = res[order(res$Start), ]
  
  return(res)
}


# Function to plot the gene tracks and labels properly:
gene.position <- function(data) {
  odd = 1
  for (i in 1:length(data$Gene)) {
    lines(x = c(data$Start[i], data$End[i]), y = c(data$Y[i], data$Y[i]), lwd = 3, col = as.character(data$Colour[i]))
    if(length(data$Gene) >= 10){
      length = abs(data$Start[i] - data$End[i])
      if(length > 5000) {
        if (odd%%2 == 0) {
          text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] - 0.1, labels = data$Gene[i], font = 3, cex = 0.7, pos = 3)
        } else {
          text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i], labels = data$Gene[i], font = 3, cex = 0.7, pos = 1)
        }
      }
    } else {
      if (odd%%2 == 0) {
        text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i] - 0.1, labels = data$Gene[i], font = 3, cex = 0.7, pos = 3)
      } else {
        text(x = (data$Start[i] + data$End[i])/2, y = data$Y[i], labels = data$Gene[i], font = 3, cex = 0.7, pos = 1)
      }
    }
    odd = odd + 1
  }
}

# Function to convert string P-value into logged P:
elog10 <- function(p) {
  if (is.character(p) & grepl('e-', p)) {
    split_p <- base::strsplit(p, split = "e-")
    tmp <- unlist(split_p)
    res <- as.numeric(tmp[2]) - log10(as.numeric(tmp[1]))
  } else {
    res <- -log10(as.numeric(p))
  }
  return(res)
}

# Function to get the LD information of specified population from the 1000 Genomes data (March 2017 release):
# NOTE: the input SNP MUST be in rsID format, not CHR:POS-based.
# NOTE: This function will leave/save the LD information in the working directory for future reference (e.g. if the user wanted to use the same LD information)
get.ld <- function(region, snp, population) {
  ld.snp = snp

  # gsub the command and filename for chr, start/end positions and the population:
  base.command = "source ~/.bashrc;
  bcftools view \
    --regions ZZ:Y1-Y2 \
    --output-type z \
    --output-file tmp.vcf.gz \
    /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/POP/POP_chrZZ.no_relatives.no_indel.biallelic.vcf.gz;

  plink1.9b4.9 \
    --vcf tmp.vcf.gz \
    --allow-no-sex \
    --snps-only \
    --r2 \
    --inter-chr \
    --ld-snp SNP \
    --ld-window-r2 0 \
    --out POP_region_ZZ.Y1-Y2_SNP;

  rm tmp.vcf.gz POP_region_ZZ.Y1-Y2_SNP.nosex"
  
  base.command = gsub(pattern = "\n ", replacement = "", base.command)
  base.command = gsub(pattern = 'ZZ', replacement = region[1], base.command)
  base.command = gsub(pattern = 'Y1', replacement = region[2], base.command)
  base.command = gsub(pattern = 'Y2', replacement = region[3], base.command)
  if(population == "TAMA"){
    base.command = gsub(pattern = 'POP/POP', replacement = "TAMAset/AFR_AMR_EAS_EUR", base.command)
    base.command = gsub(pattern = 'POP', replacement = population, base.command)
  } else{
    base.command = gsub(pattern = 'POP', replacement = population, base.command)
  }
  base.command = gsub(pattern = 'SNP', replacement = ld.snp, base.command)
  
  # Make a system call to run the bcftools/plink command.
  # I'm only assigning it to a variable to suppress any possible form of
  # messages/outputs from the command, just in case
  # ignores any errors when running the LD command, but will output the error to your screen
  messages = system(base.command, ignore.stdout = TRUE, intern = TRUE, ignore.stderr = TRUE)
  
  # Import the LD data:
  ld.file = "POP_region_ZZ.Y1-Y2_SNP.ld"
  ld.file = gsub(pattern = 'ZZ', replacement = region[1], ld.file)
  ld.file = gsub(pattern = 'Y1', replacement = region[2], ld.file)
  ld.file = gsub(pattern = 'Y2', replacement = region[3], ld.file)
  ld.file = gsub(pattern = 'POP', replacement = population, ld.file)
  ld.file = gsub(pattern = 'SNP', replacement = ld.snp, ld.file)
  
  # Check the ld file was made & import it
  if(ld.file %in% list.files(pattern = ".ld")){
    ld = read.table(ld.file, stringsAsFactors = FALSE, header = TRUE)
  } else{
    message("Top SNP / specified SNP not in 1000 Genomes biallelic SNPs")
    ld = data.frame(CHR_A = NA, BP_A = NA, SNP_A = NA, CHR_B = NA, BP_B = NA, SNP_B = NA, R2 = NA)
  }
  # return the ld file data
  return(ld)
}

# Function to read in and pull out relevant info from PLINK clump output:
read.plink.loci <- function(file = NULL) {
  if (is.null(file)) {
    stop('You must provide a file for reading')
  }
  data = read.table(file, stringsAsFactors = FALSE, header = TRUE)
  data = data[,c(1,3:5)]
  return(data)
}
