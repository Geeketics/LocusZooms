######################
# Locus Zoom, Make LD Data
# July 2020
# Tanya Major & Riku Takei
# Uni of Otago

#### Important Running Notes: ####

### Compulsory flags:

## One of snp, gene, or region must be specified to create the plot:
# snp: specify the SNP to be annotated (you must also include ignore.lead = TRUE if choosing this option)
# gene: specify the Gene to make the plot around
# region: specify the chromsome region you want to plot (must be specified as c(chr, start, end)

## As well as each of the following:
# data: specify the data.frame (or a list of data.frames) to be used in the plot (requires the columns "CHR", "BP", "SNP", and either "P" or "logBF")
# genes.data: specify a data.frame with gene locations to plot beneath the graph (requires the columns "Gene", "Chrom", "Start", "End", and "Coding") # the Gencode or UCSC {Gencode,UCSC}_GRCh37_Genes_UniqueList{2017,2021}.txt files in this repo can be used for this
# plot.title: specify a title to go above your plot
# file.name: specify a filename for your plot to be saved to

### Optional flags:

# ld.file: specify a data.frame with LD values relevant to the SNP specified by snp (requires the columns "SNP_B" and "R2") 
# offset_bp: specify how far either side of the snp, gene, or region you want the plot to extend (defaults to 200000)
# psuedogenes: when using one of the three gene lists in this repo you can specify whether you want to plot the pseudogenes (defaults to FALSE)
# RNAs: when using one of the two gene lists created in 2021 in this repo you can specify whether you want to plot lncRNA and ncRNA genes (defaults to FALSE)
# plot.type: specify the file format of the plot (defaults to "jpg", options are "jpg", "svg", or "view_only" which will not save the plot, but output it to RStudio Viewer instead)
# nominal: specify the nominal significance level to draw on the plot (in -log10(_P_), default is 6 or _P_ = 1e-6)
# significant: specify the significance level to draw on the plot (in -log10(_P_), default is 7.3 or _P_ = 5e-8) 
# secondary.snp: provide the list of secondary SNP IDs (must match IDs in results file) to be highlighted on the plot
# secondary.label: specify whether to label the secondary SNPs on the plot (defaults to FALSE)
# secondary.circle: specify whether to add a red circle around the secondary SNPs on the plot (defaults to TRUE)
# genes.pvalue: specify a data.frame of p-values (e.g. MAGMA results) associated with each gene (requires the columns "Gene" and "P") 
# colour.genes: specify whether to colour genes based on a p-value provided in gene.pvalue (defaults to FALSE)
# population: specify the 1000 genomes population to use when calculating LD if ld.file = NULL (defaults to "EUR", options are "AFR", "AMR", "EAS", "EUR", "SAS", "TAMA", and "ALL")
# sig.type: specify whether the y-axis should be labelled as -log10(_P_) or log10(BF) (defaults to "P", options are "P", "logP", or "logBF"). For the "P" option an additional -log10 conversion of the input "P" column will be performed.
# nplots: specify whether multiple results plots will be saved into your jpeg file (e.g. plot two GWAS results one above another; defaults to FALSE)
# ignore.lead: specify whether to ignore the SNP with the smallest P and use the SNP specified by 'snp' to centre the plot (defaults to FALSE)
# rsid.check: specify whether to check if the SNPs are labelled with rsIDs # should only matter if script is calculating LD for you (defaults to TRUE)
# nonhuman: specify whether the data to plot has come from a non-human sample-set (defaults to FALSE) # if the data going in is from a non-human species make sure the chromosome column is only numbers (e.g. 1 instead of chr1, 23 instead of X). 
                                                           


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

#### Function to make LocusZoom like plots ####
locus.zoom <- function(data = NULL, snp = NA, gene = NA, region = NA, ld.file = NULL, offset_bp = 200000, genes.data = NULL, psuedogenes = FALSE, RNAs = FALSE, plot.title = NULL, plot.type = "jpg", nominal = 6, significant = 7.3, file.name = NULL, secondary.snp = NA, secondary.label = FALSE, secondary.circle = TRUE, genes.pvalue = NULL, colour.genes = FALSE, population = "EUR", sig.type = "P", nplots = FALSE, ignore.lead = FALSE, rsid.check = TRUE, nonhuman = FALSE) {
  
  # Define constants:
  LD.colours <- data.frame(LD = as.character(seq(from = 0, to = 1, by = 0.1)), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2)), stringsAsFactors = FALSE)
  GENE.colours <- data.frame(Threshold = c(">2.6e-6", ">1e-10", ">1e-15", ">1e-20", "<1e-20"), Colour = c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), stringsAsFactors = FALSE)
  
  # load scales library
  if("scales" %in% data.frame(installed.packages())[, "Package"]){
    library("scales")
  } else{
    stop("This function requires the package 'scales' to run.\nUse install.packages('scales') to install the package before running this code again.")
  }
  
  
  # Load Data
  # If plotting multiple summary stats, take the first summary stats as lead/reference data:
  if (is.data.frame(data)) {
    lead.data = data
  } else {
    lead.data = data[[1]]
  }
  
  # Error check data header
  if(sig.type == "P") {
    if(!all(c("CHR", "BP", "SNP", "P") %in% names(lead.data))){
      stop("Your data file does not contain a CHR, BP, SNP, or P column.\nCheck your header line.")
    }
  } else if(sig.type == "logP") {
    if(!all(c("CHR", "BP", "SNP", "logP") %in% names(lead.data))){
      stop("Your data file does not contain a CHR, BP, SNP, or logP column.\nCheck your header line.")
    }
  } else if(sig.type == "logBF") {
    if(!all(c("CHR", "BP", "SNP", "logBF") %in% names(lead.data))){
      stop("Your data file does not contain a CHR, BP, SNP, or logBF column.\nCheck your header line.")
    }
  } else {
      stop("Unrecognised significance type. The options are P, logP, or logBF.")
  }

  lead.data$SNP = as.character(lead.data$SNP)
  
  # Check the SNPs are in rsID format and no duplicates:
  if (rsid.check) {
    check.rsid(lead.data$SNP)
  }
  
  # Check gene data header
  if(!all(c("Chrom", "Start", "End", "Coding") %in% names(genes.data))){
    stop("Your genes.data file does not contain a Chrom, Start, End, or Coding column.\nCheck your header line.")
  }
  
  # convert 'chr1' to '1' 
  if(!nonhuman){
    if(!(class(genes.data$Chrom) %in% c("integer", "numeric"))){
      genes.data$Chrom <- gsub(genes.data$Chrom, pattern = "chr", replacement = "")
      genes.data$Chrom[genes.data$Chrom == "X"] <- 23
      genes.data$Chrom[genes.data$Chrom == "Y"] <- 24
      genes.data$Chrom[genes.data$Chrom %in% c("M", "MT")] <- 26
      genes.data$Chrom <- as.numeric(genes.data$Chrom)
    }
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
  
  ## Pull out the relevant information from the gene data.

  # Any gene that overlaps/intersect with the defined region is included:
  genes.data = genes.data[genes.data$Chrom == region[1], ]
  genes.data = genes.data[genes.data$End > region[2], ]
  genes.data = genes.data[genes.data$Start < region[3], ]
  
  # Remove psuedogenes & RNA Info:
  if(!psuedogenes) {
    genes.data = genes.data[!(genes.data$Coding %in% c("psuedogene", "Non-Coding")), ]
  }
  
  if(!RNAs) {
    genes.data = genes.data[!(genes.data$Coding %in% c("lncRNA", "ncRNA")), ]
  }
  
  # Pull out the relevant information from the gene p-values data
  if(colour.genes) {
    if(!all(c("Gene", "P") %in% names(genes.pvalue))){
      stop("Your genes.pvalue file does not contain a Gene or P column.\nCheck your file header.")
    }
    genes.pvalue = genes.pvalue[genes.pvalue$Gene %in% genes.data$Gene, ]
  }
    
  # Pull out the relevant data from the result file(s), and logBF/log p-value:
  if (is.data.frame(data)) {
    data = subset.data(data, region)
    if (sig.type == "P") {
      data$logP = as.numeric(unlist(lapply(data$P, elog10)))
    } else if (sig.type == "logBF") {
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
    if(length(lead.ind) == 0){
      stop(paste0("Your SNP (", snp, ") is not present in your data."))
    }
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
  
### Make Plot ###
  
  # Define output plot size
  npanel = ifelse(nplots, length(data), 1)
  plot.height = (npanel * 80) + 50
  if(plot.type == "jpg"){
    jpeg(width = 160, height = plot.height, units = "mm", res = 300, filename = file.name)
  } else if(plot.type == "svg") {
    svg(width = (160 / 25.4), height = (plot.height / 25.4), filename = file.name)
  } else if(plot.type != "view_only"){
    stop("Unrecognised plot type. The options are jpg, png, or view_only.")
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
      plot.locus(data.plot = tmp.dat, plot.title = names(data.plot)[i], secondary.snp = secondary.snp, secondary.label = secondary.label, secondary.circle = secondary.circle, sig.type = sig.type, plot.var = plot.var)
      rm(tmp.dat)
    }
  } else {
    y.max = max(data.plot$logP, 8)
    plot.var = c(y.max, x.min, x.max, lead.snp, nominal, significant)
    plot.locus(data.plot = data.plot, plot.title = plot.title, secondary.snp = secondary.snp, secondary.label = secondary.label, secondary.circle = secondary.circle, sig.type = sig.type, plot.var = plot.var)
  }
  
  # Plot Gene tracks
  par(mar = c(4, 4, 0.5, 8), mgp = c(2, 1, 0), xpd = FALSE)
  if(length(genes.data[, "Gene"]) > 15){
    track.max = 6
    font.size = 0.45
  } else{
    track.max = 3
    font.size = 0.6
  }
  plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", lead.chr), ylab="", xlim = c(x.min, x.max), ylim = c(0, track.max), xaxt = "n")
  x_marks = axTicks(side = 1)
  axis(side = 1, at = x_marks, labels = format(x_marks, scientific = FALSE, big.mark = ",", trim = TRUE))

  if (nrow(genes.data) != 0) {
                       
  # add colour column to genes.data
  if(colour.genes) {
    genes.data = merge.gene.colour(genes.data, genes.pvalue, GENE.colours)
  } else {
    genes.data$Colour = "#7F7F7F"
  }
  
  # Stagger the genes
    if(length(genes.data[, "Gene"]) > 15){
      y = rep(c(5, 2, 3, 4, 1), times = length(genes.data[ ,"Gene"]))
      genes.data$Y = y[1:length(genes.data$Gene)]
    } else{
      y = rep(c(2.5, 1.5, 0.5), times = length(genes.data[ ,"Gene"]))
      genes.data$Y = y[1:length(genes.data$Gene)]
    }
    # Plot the gene tracks:
    for(track in unique(genes.data$Y)){
      genes.set <- genes.data[genes.data$Y == track, ]
      if (nrow(genes.set) > 0) {
        gene.position(genes.set, fontsize = font.size, plot.var = plot.var)
      }
    }
  }

  # add gene colour legend
  if(colour.genes) {
    legend.colour = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080", "#7F7F7F")
    par(xpd = TRUE)
    legend(x = "right", legend = c(expression("<1x10"^-20), expression(paste("<1x10"^-15, "; ≥1x10"^-20)), expression(paste("<1x10"^-10, "; ≥1x10"^-15)), expression(paste("<2.6x10"^-6, "; ≥1x10"^-10)), expression("≥2.6x10"^-6), "Unknown"), col = legend.colour, fill = legend.colour, border = legend.colour, pt.cex = 1.2, cex = 0.8, bg = "white", box.lwd = 0, title = "p-value", inset = -0.22)
  }

  if( plot.type != "view_only" ) {
    dev.off()
  }
}

#### Function to create LocusZoom style plot (without gene track): ####
plot.locus <- function(data.plot = NULL, plot.title = NULL, secondary.snp = NA, secondary.label = FALSE, secondary.circle = TRUE, sig.type = "P", plot.var = NULL) {
  # Variables:
  y.max = as.numeric(plot.var[1])
  x.min = as.numeric(plot.var[2])
  x.max = as.numeric(plot.var[3])
  lead.snp = plot.var[4]
  nominal = as.numeric(plot.var[5])
  significant = as.numeric(plot.var[6])
  
  # Plot SNP presence:
  par(mar = c(0, 4, 2, 8), mgp = c(2, 1, 0), xpd = FALSE)
  plot(x = data.plot$BP, y = rep(1, times = nrow(data.plot)), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted\nSNPs", las = 2, xlim = c(x.min, x.max), cex.lab = 0.8, col = alpha(colour = "black", alpha = 0.2))
  title(plot.title, line = 0)
  
  # Plot Manhattan/LocusZoom of region
  par(mar = c(0, 4, 0, 8), mgp = c(2, 1, 0), xpd = FALSE)
  ylab = ifelse(sig.type == "P" | sig.type == "logP", expression(-log[10](italic(P))), expression(log[10](BF)))
  plot(x = data.plot$BP, y = data.plot$logP, ylim = c(0, y.max*1.1), pch = 20, col = as.character(data.plot$Colour), xlab = "", ylab = ylab, cex = 0.8, xaxt = "n", xlim = c(x.min, x.max))
  abline(h = nominal, col = "blue", lty = "dashed")
  abline(h = significant, col = "red", lty = "dashed")
  
  # Plot the lead SNP
  if (lead.snp %in% data.plot$SNP) {
    ind = which(data.plot$SNP == lead.snp)
    lead.pos = data.plot$BP[ind]
    lead.logp = data.plot$logP[ind]
    points(x = lead.pos, y = lead.logp, pch = 18, cex = 1.5, col = "#7D26CD")
    text(x = lead.pos, y = lead.logp, labels = lead.snp, pos = 3)
  }
  
  # Plot label/text for the secondary SNP (removes lead.snp from the list of secondary SNPs before plotting)
  if(any(!is.na(secondary.snp))){
    secondary.snp <- secondary.snp[secondary.snp != lead.snp]
    check = which(data.plot$SNP %in% secondary.snp)
    if (length(check) != 0) {
      secondary.data = data.plot[data.plot$SNP %in% c(secondary.snp, lead.snp), ]
      plot.secondary.point(data = secondary.data, snps = secondary.data$SNP, lead.snp = lead.snp, plot.data = data.plot, plot.var = plot.var, label = secondary.label, circle = secondary.circle)
    }
  }
  
  # Add LD legend
  legend.colour = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080", "#7F7F7F")
  par(xpd = TRUE)
  legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2", "Unknown"), col = legend.colour, fill = legend.colour, border = legend.colour, pt.cex = 1.2, cex = 0.8, bg = "white", box.lwd = 0, title = expression("r"^2), inset = c(-0.14, 0.01))
}

##### Function to check if the input variants have rsIDs ####
check.rsid <- function(snp = NULL) {
  # Stop if CHR:POS ID:
  if (any(!grepl('rs', snp))) {
    stop("Your SNP column does not have rsIDs")
  }
  # Stop if there are duplicate SNPs:
  if (length(which(duplicated(snp))) > 0) {
    stop("There are duplicate rsIDs in your results file - Please remove them before running again")
  }
}

#### Function to make regions out of variant or gene information ####
get.region <- function(snp.dat, snp, gene.dat, gene) {
  # If SNP is given:
  if (!is.na(snp)) {
    snp.ind = which(snp.dat$SNP == snp)
    
    if(length(snp.ind) == 0){
      stop(paste0("Your SNP (", snp, ") is not present in your data."))
    }
    
    snp.chr = snp.dat$CHR[snp.ind]
    snp.pos = snp.dat$BP[snp.ind]
    region = c(snp.chr, snp.pos, snp.pos)
  } else if (!is.na(gene)) {
  # If Gene is given
    gene.ind = which(gene.dat$Gene == gene)
    
    if(length(gene.ind) == 0){
      stop(paste0("Your gene (", gene, ") is not present in your data."))
    }
    
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

#### Function to subset relevant region from a dataset ####
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

#### Function to merge the LD and LD colours with the relevant region of the summary stats: ####
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

#### Function to plot secondary SNPs: ####
plot.secondary.point <- function(data, snps, lead.snp, plot.data, plot.var, label = FALSE, circle = TRUE) {
  
  lead.ind = which(data$SNP == lead.snp)
  lead.pos = data$BP[lead.ind]
  
  data = data[which(data$SNP != lead.snp), ]
  data = data[data$BP >= as.numeric(plot.var[2]) & data$BP <= as.numeric(plot.var[3]), ]
  
  
  # plot red line around secondary SNP
  if (circle) {
    points(x = data$BP, y = data$logP, pch = 1, cex = 1.1, col = "#FF0000")
  }
  
  # add SNP labels if requested
  if (label) {
    # set up labeling offsets (x-axis)
    x.min = as.numeric(plot.var[2])
    x.max = as.numeric(plot.var[3])
    x.offset = abs(x.max - x.min) / 150 * 15
    
    data$label.x.offset[data$BP < lead.pos] = data$BP[data$BP < lead.pos] - x.offset / 3
    data$side[data$BP < lead.pos] = 2
    
    data$label.x.offset[data$BP >= lead.pos] = data$BP[data$BP >= lead.pos] + x.offset / 3
    data$side[data$BP >= lead.pos] = 4
    
    data$label.y.offset = NA
    for(snp in data$SNP){
      ind = which(data$SNP == snp)
      pos = data$BP[ind]
      
      # set up labeling offsets (y-axis) - considers SNPs around it  
      surrounding.data = plot.data[plot.data$BP > (pos - (abs(x.max - x.min) / 6)) & plot.data$BP < (pos + (abs(x.max - x.min) / 6)), ]
      surrounding.data = surrounding.data[surrounding.data$logP > data$logP[ind] & surrounding.data$logP < data$logP[ind] + 5, ]
      if(length(surrounding.data[, 1]) == 0){
        data$label.y.offset[ind] = (data$logP[ind] + 2) * 1.03
      } else {
        data$label.y.offset[ind] = max(surrounding.data$logP) * 1.03
      }
      rm(surrounding.data, ind, pos)
    }
    
    for(offset.pos in unique(data$label.y.offset)){
      if(length(data[data$label.y.offset == offset.pos, 1]) > 1){
        data$label.y.offset[data$label.y.offset == offset.pos] = jitter(data$label.y.offset[data$label.y.offset == offset.pos], amount = 1)
      }
    }
    
    for(snp in data$SNP){
      ind = which(data$SNP == snp)
      pos = data$BP[ind]
      
      logp = data$logP[ind]
      label.x = data$label.x.offset[ind]
      label.y = data$label.y.offset[ind]
      side = data$side[ind]
      
      # add lines and text to label SNP
      text(x = label.x, y = (label.y * 1.01), labels = snp, cex = 0.7, pos = side, offset = 0.2)
      segments(x0 = pos, x1 = label.x, y0 = logp, y1 = label.y * 1.01)
    }
  }
}

#### Function to merge the gene and gene colours with the relevant region of the summary stats: ####
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


#### Function to plot the gene tracks and labels properly: ####
gene.position <- function(data, fontsize = 0.6, plot.var = NULL) {
  # Variables:
  x.min = as.numeric(plot.var[2])
  x.max = as.numeric(plot.var[3])
  plot.length = x.max - x.min
  
  # narrow down to genes to be labeled before working out whether to label above/below line
  if(plot.length > 5000000) {
    gene.length = 16000
  } else if(plot.length > 2000000) {
    gene.length = 12000
  } else {
    gene.length = 8000
  }

  if(length(data$Gene) >= 10) {
    text_data <- data[abs(data$Start - data$End) > gene.length, ]
  } else{
    text_data <- data
  }
  
  # add gene lines & labels to plot
  odd = 1
  for (i in 1:length(data$Gene)) {
    lines(x = c(data$Start[i], data$End[i]), y = c(data$Y[i], data$Y[i]), lwd = 3, col = as.character(data$Colour[i]))
  }
  
  for(i in 1:length(text_data$Gene)) {
    mid.point = (max(x.min, text_data$Start[i]) + min(text_data$End[i], x.max))/2
    if (odd%%2 == 0) {
      text(x = mid.point, y = text_data$Y[i], labels = text_data$Gene[i], font = 3, cex = fontsize, pos = 3, offset = 0.2)
    } else {
      text(x = mid.point, y = text_data$Y[i], labels = text_data$Gene[i], font = 3, cex = fontsize, pos = 1, offset = 0.25)
    }
    odd = odd + 1
  }
}


#### Function to convert string P-value into logged P: ####
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

#### Function to get the LD information of specified population from the 1000 Genomes data (March 2017 release): ####
# NOTE: the input SNP MUST be in rsID format, not CHR:POS-based.
# NOTE: This function will leave/save the LD information in the working directory for future reference (e.g. if the user wanted to use the same LD information)
get.ld <- function(region, snp, population) {
  ld.snp = snp

  if (region[1] == "23") {
    region[1] = "X"
  }
  
  vcf.filename = "POP_chrZZ.no_relatives.no_indel.biallelic.vcf.gz"
  vcf.filename = gsub(pattern = 'ZZ', replacement = region[1], vcf.filename)
  if(population == "TAMA"){
    vcf.filename = gsub(pattern = 'POP', replacement = "AFR_AMR_EAS_EUR", vcf.filename)
  } else{
    vcf.filename = gsub(pattern = 'POP', replacement = population, vcf.filename)
  }
  
  # check necessary 1000 genomes file can be reached
  if(!(vcf.filename %in% list.files(path = paste0("/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/", population)))){
    stop(paste0("The file ", paste0("/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/", population, "/", vcf.filename), " cannot be found."))
  }
  
  # gsub the command and filename for chr, start/end positions and the population:
  base.command = "source ~/.bashrc;
  bcftools view \
    --regions ZZ:Y1-Y2 \
    --output-type z \
    --output-file tmp.vcf.gz \
    /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/POP/1000VCF;

  plink \
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
  base.command = gsub(pattern = "1000VCF", replacement = vcf.filename, base.command)
  base.command = gsub(pattern = 'POP', replacement = population, base.command)
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

#### Function to read in and pull out relevant info from PLINK clump output: ####
read.plink.loci <- function(file = NULL) {
  if (is.null(file)) {
    stop('You must provide a file for reading')
  }
  data = read.table(file, stringsAsFactors = FALSE, header = TRUE)
  data = data[,c(1,3:5)]
  return(data)
}
