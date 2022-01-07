######################
# Locus Zoom, Make LD Data
# July 2020
# Tanya Major & Riku Takei
# Uni of Otago

### Important Running Notes:
## data: expects a data.frame (or a list of data.frames) with at least columns containing the chromosome, positions, rsIDs, and p-values of your results - these must be labelled CHR, BP, SNP, and P (or logBF)
## snp: specifies a SNP to centre the plot around
## gene: specifies a gene to centre the plot around
## region: specifies the chromosome start and end you wish to plot
## ld.file: expects a data.frame of the LD between your lead SNP and all other SNPs (requires the columns SNP_B and R2) - if left blank the script can calculate this for you so long as you have access to the biochem servers
## offset_bp: specifies how far either side of your gene/snp/region of interest to plot (in base pairs)
## genes.data: specifies a data.frame of genes within the plot region - expects the headers Gene, Chrom, Start, End, & Coding
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
## sig.type: specifies whether the y-axis should be -log10(P) or -log10(BF) - the options are P (will be converted to -log10(P)), logP (no conversion needed), or logBF.
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
