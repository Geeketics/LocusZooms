######################
# Locus Zoom, Make LD Data
# September 2016/February 2017
# Tanya Flynn
# Uni of Otago

### Important Running Notes:
## You need to have acces to the data that the associations were originally done on
## Genes.Data expects a data.frame with headers Gene, Chrom, Start, & End - this will be used to annotate the graph.

## LD Colours
# 1.0 - 0.8 = #FF0000
# 0.8 - 0.6 = #FFA500
# 0.6 - 0.4 = #00FF00
# 0.4 - 0.2 = #87CEFA
# 0.2 - 0.0 = #000080
# NA = #7F7F7F
# top-hit = #7D26CD

locus.zoom <- function(CHR = NULL, BP = NULL, P = NULL, SNP.List = NULL, SNP = NA, Gene = NA, LD.File = NULL, kb = NULL, Genes.Data = NULL, NonCoding = FALSE, Plot.Title = NULL, Nominal = 6, Significant = 7.3, File.Name = NULL){
  # Load Data
  results.data <- data.frame(snps = SNP.List, chr = CHR, pos = BP, p = P)
  results.data[,"snps"] <- as.character(results.data[,"snps"])
  results.data <- results.data[!duplicated(results.data$snps),]
  ld.data <- LD.File
  genes.data <- Genes.Data

  # Reduce data to relevant segment
  if(exists("kb") == FALSE){
    kb <- 200000
  }
  
  if(is.na(SNP) & is.na(Gene)){
    stop("You must specify a SNP or a Gene to plot")
  }
  
  if(!is.na(Gene)){
    gene.start <- genes.data[genes.data$Gene == Gene, "Start"]
    gene.end <- genes.data[genes.data$Gene == Gene, "End"]
    Chr <- genes.data[genes.data$Gene == Gene, "Chrom"]
    new.results.data <- results.data[results.data$chr == Chr & results.data$pos >= (gene.start - kb) & results.data$pos <= (gene.end + kb),]
    p.min <- min(new.results.data[,"p"], na.rm = TRUE)
    SNP <- new.results.data[new.results.data[,"p"] == p.min & !is.na(new.results.data[,"p"]), "snps"]
    snp.pos <- new.results.data[new.results.data[,"snps"] == SNP, "pos"]
    snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SNP, "p"])
    genes.data <- genes.data[genes.data$Chrom == Chr & genes.data$End > (gene.start - kb) & genes.data$Start < (gene.end + kb),]
    y = rep(c(1.5, 0.5), times = length(genes.data[,"Gene"]))
    genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
  } else{
    snp.pos <- results.data[results.data[,"snps"] == SNP, "pos"]
    Chr <- results.data[results.data[,"snps"] == SNP, "chr"]
    new.results.data <- results.data[results.data$chr == Chr & results.data$pos >= (snp.pos - kb) & results.data$pos <= (snp.pos + kb),]
    p.min <- min(new.results.data[,"p"], na.rm = TRUE)
    snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SNP, "p"])
    genes.data <- genes.data[genes.data$Chrom == Chr & genes.data$End > (snp.pos - kb) & genes.data$Start < (snp.pos + kb),]
    y = rep(c(1.5, 0.5), times = length(genes.data[,"Gene"]))
    genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
  }

  # Remove Non-Coding Gene Info
  if(NonCoding == FALSE){
    genes.data <- genes.data[genes.data$Coding != "Non-Coding",]
  }
  
  # Add LD to Results
  new.ld.data <- ld.data[,c("snps", SNP)]

  round.up <- function(x, decimals = 1){
    round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
  }

  new.results.data <- merge(new.results.data, new.ld.data, by = "snps", all.x = TRUE)
  new.results.data[,"plot.ld"] <- round.up(new.results.data[,SNP], decimals = 1)
  new.results.data[new.results.data[,"plot.ld"] > 1 & !is.na(new.results.data[,"plot.ld"]), "plot.ld"] = 1
  
  LD.colours <- data.frame(LD = c(seq(from = 0, to = 1, by = 0.1), "NaN"), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2), "#7F7F7F"))
  
  new.results.data.plot <- merge(new.results.data, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)
  genes.top <- genes.data[genes.data$Y == 1.5,]
  genes.bot <- genes.data[genes.data$Y == 0.5,]
  
  # Make Plotting Variables
  p.max <- -log10(p.min)
  y.max <-  max(round.up(p.max, decimals = 0), 8)
  x.min <-  (snp.pos - kb)
  x.max <-  (snp.pos + kb)
  text_offset = abs(x.max - x.min)/150 * 15
  
  # Make Plot
  jpeg(width = 150, height = 150, units = "mm", res = 300, file = File.Name)
  layout(matrix(c(1, 2, 3), byrow = TRUE), heights = c(2, 6, 3))
  
  ## plot 1 - where SNPs are
  par(mar = c(0.6, 4, 4, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = rep(1, times = length(new.results.data.plot[,"snps"])), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted SNPs", las = 2, main = Plot.Title, xlim = c(x.min, x.max))
  
  ## plot 2 - actual manhattan p-value
  par(mar = c(0.5, 4, 0.6, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = -log10(new.results.data.plot[,"p"]), ylim = c(0, y.max), pch = 20, col = as.character(new.results.data.plot[,"Colour"]), xlab = "", ylab = expression(-log[10](italic(P))), cex = 1.5, xaxt = "n", xlim = c(x.min, x.max))
  points(x = snp.pos, y = snp.p, pch = 18, cex = 2, col = "#7D26CD")
  text(x = (snp.pos + text_offset), y = p.max, labels = SNP)
  abline(h = Nominal, col = "blue", lty = "dashed")
  abline(h = Significant, col = "red", lty = "dashed")
  legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2"), col = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), fill = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), border = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), pt.cex = 2, cex = 1.2, bg = "white", box.lwd = 0, title = expression("r"^2))
  
  ## plot 3 - where the genes are
  par(mar = c(4, 4, 0.5, 4), mgp = c(2, 1, 0))
  plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", Chr), ylab="", xlim = c(x.min, x.max), ylim = c(0,2))
  odd = 1
  for(gene in 1:length(genes.top[,"Gene"])){
    lines(x = c(genes.top[gene,"Start"], genes.top[gene,"End"]), y = c(genes.top[gene,"Y"], genes.top[gene,"Y"]), lwd = 3, col = "#000080")
    if(length(genes.data[,"Gene"]) > 10){
      length = abs(genes.top[gene,"Start"] - genes.top[gene,"End"])
      if(length > 8000 & odd%%2 == 0){ ## using the %% tells you if a number is odd/even odd == 1, even == 0
        text(x = (genes.top[gene,"Start"] + genes.top[gene,"End"])/2, y = genes.top[gene, "Y"] - 0.2, labels = genes.top[gene, "Gene"], cex = 0.8)
      } else{
        if(length > 8000 & odd%%2 == 1){
          text(x = (genes.top[gene,"Start"] + genes.top[gene,"End"])/2, y = genes.top[gene, "Y"] + 0.2, labels = genes.top[gene, "Gene"], cex = 0.8)
        }
      }
    } else{
      text(x = (genes.top[gene,"Start"] + genes.top[gene,"End"])/2, y = genes.top[gene, "Y"] + 0.2, labels = genes.top[gene, "Gene"])
    }
      odd = odd + 1
  }
  odd = 1
  for(gene in 1:length(genes.bot[,"Gene"])){
    lines(x = c(genes.bot[gene,"Start"], genes.bot[gene,"End"]), y = c(genes.bot[gene,"Y"], genes.bot[gene,"Y"]), lwd = 3, col = "#000080")
    if(length(genes.data[,"Gene"]) > 10){
      length = abs(genes.bot[gene,"Start"] - genes.bot[gene,"End"])
      if(length > 8000 & odd%%2 == 0){ ## using the %% tells you if a number is odd/even odd == 1, even == 0
        text(x = (genes.bot[gene,"Start"] + genes.bot[gene,"End"])/2, y = genes.bot[gene, "Y"] - 0.2, labels = genes.bot[gene, "Gene"], cex = 0.8)
      } else{
        if(length > 8000 & odd%%2 == 1){
          text(x = (genes.bot[gene,"Start"] + genes.bot[gene,"End"])/2, y = genes.bot[gene, "Y"] + 0.2, labels = genes.bot[gene, "Gene"], cex = 0.8)
        }
      }
    } else{
      text(x = (genes.bot[gene,"Start"] + genes.bot[gene,"End"])/2, y = genes.bot[gene, "Y"] + 0.2, labels = genes.bot[gene, "Gene"])
    }
    odd = odd + 1
  }
  dev.off()

rm(Genes.Data, gene, genes.data, LD.colours, ld.data, LD.File, new.ld.data, new.results.data, new.results.data.plot, results.data, BP, Chr, CHR, kb, Nominal, P, p.max, p.min, Plot.Title, Significant, SNP, SNP.List, snp.p, snp.pos, x.max, x.min, y, y.max, text_offset)
  
}
