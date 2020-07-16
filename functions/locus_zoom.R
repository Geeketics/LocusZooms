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

locus.zoom <- function(CHR = NULL, BP = NULL, P = NULL, SNP.List = NULL, SNP = NA, Gene = NA, Region = NA, LD.File = NULL, basepairs = 200000, Genes.Data = NULL, NonCoding = FALSE, Plot.Title = NULL, Nominal = 6, Significant = 7.3, File.Name = NULL, SecondarySNP = NA){
  # region specified as c(chr, start, end)
  # Load Data
  results.data <- data.frame(snps = SNP.List, chr = CHR, pos = BP, p = P)
  results.data[,"snps"] <- as.character(results.data[,"snps"])
  results.data <- results.data[!duplicated(results.data$snps),]
  ld.data <- LD.File
  genes.data <- Genes.Data

  # Reduce data to relevant segment
  if(!is.na(Gene)){
    gene.start <- genes.data[genes.data$Gene == Gene, "Start"]
    gene.end <- genes.data[genes.data$Gene == Gene, "End"]
    Chr <- genes.data[genes.data$Gene == Gene, "Chrom"]
    new.results.data <- results.data[results.data$chr == Chr & results.data$pos >= (gene.start - basepairs) & results.data$pos <= (gene.end + basepairs),]
    p.min <- min(new.results.data[,"p"], na.rm = TRUE)
    SNP <- new.results.data[new.results.data[,"p"] == p.min & !is.na(new.results.data[,"p"]), "snps"]
    snp.pos <- new.results.data[new.results.data[,"snps"] == SNP, "pos"]
    snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SNP, "p"])
    genes.data <- genes.data[genes.data$Chrom == Chr & genes.data$End > (gene.start - basepairs) & genes.data$Start < (gene.end + basepairs),]
    y = rep(c(1.5, 0.5), times = length(genes.data[,"Gene"]))
    genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
  } else{
    if(!is.na(SNP)){
      snp.pos <- results.data[results.data[,"snps"] == SNP, "pos"]
      Chr <- results.data[results.data[,"snps"] == SNP, "chr"]
      new.results.data <- results.data[results.data$chr == Chr & results.data$pos >= (snp.pos - basepairs) & results.data$pos <= (snp.pos + basepairs),]
      p.min <- min(new.results.data[,"p"], na.rm = TRUE)
      snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SNP, "p"])
      genes.data <- genes.data[genes.data$Chrom == Chr & genes.data$End > (snp.pos - basepairs) & genes.data$Start < (snp.pos + basepairs),]
    } else{
      if(!is.na(Region)[1]){
        Chr <- Region[1]
        gene.start <- Region[2]
        gene.end <- Region[3]
        new.results.data <- results.data[results.data$chr == Chr & results.data$pos >= (gene.start - basepairs) & results.data$pos <= (gene.end + basepairs),]
        p.min <- min(new.results.data[,"p"], na.rm = TRUE)
        SNP <- new.results.data[new.results.data[,"p"] == p.min & !is.na(new.results.data[,"p"]), "snps"]
        snp.pos <- new.results.data[new.results.data[,"snps"] == SNP, "pos"]
        snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SNP, "p"])
        genes.data <- genes.data[genes.data$Chrom == Chr & genes.data$End > (gene.start - basepairs) & genes.data$Start < (gene.end + basepairs),]
        y = rep(c(1.5, 0.5), times = length(genes.data[,"Gene"]))
        genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
      } else{
        stop("You must specify a SNP, Gene or Region to plot")
      }
    }
  }

  # Remove Non-Coding Gene Info
  if(NonCoding == FALSE){
    genes.data <- genes.data[genes.data$Coding != "Non-Coding",]
  }
  
  # Add LD to Results
  new.ld.data <- ld.data[,c("SNP_B", "R2")]

  round.up <- function(x, decimals = 1){
    round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
  }

  new.results.data <- merge(new.results.data, new.ld.data, by.x = "snps", by.y = "SNP_B", all.x = TRUE)
  new.results.data[, "R2"] <- as.numeric(new.results.data[, "R2"])
  new.results.data[,"plot.ld"] <- round.up(new.results.data[,"R2"], decimals = 1)
  new.results.data[new.results.data[,"plot.ld"] > 1 & !is.na(new.results.data[,"plot.ld"]), "plot.ld"] = 1
  new.results.data[, "plot.ld"] <- as.character(new.results.data[, "plot.ld"])
  
  LD.colours <- data.frame(LD = as.character(seq(from = 0, to = 1, by = 0.1)), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2)), stringsAsFactors = FALSE)
  
  new.results.data.plot <- merge(new.results.data, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)
  # make all variants without an LD value grey
  new.results.data.plot[is.na(new.results.data.plot[, "Colour"]), "Colour"] <- "#7F7F7F"
  # sort file based on position
  new.results.data.plot <- new.results.data.plot[order(new.results.data.plot$pos), ]
  plot.last <- new.results.data.plot[new.results.data.plot[, "plot.ld"] > 0.2 & !is.na(new.results.data.plot[, "plot.ld"]), ]
  new.results.data.plot <- new.results.data.plot[new.results.data.plot[, "plot.ld"] <= 0.2 | is.na(new.results.data.plot[, "plot.ld"]), ]
  new.results.data.plot <- rbind(new.results.data.plot, plot.last)
  rm(plot.last)

  y = rep(c(2.5, 1.5, 0.5), times = length(genes.data[,"Gene"]))
  genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
  genes.top <- genes.data[genes.data$Y == 2.5,]
  genes.mid <- genes.data[genes.data$Y == 1.5,]
  genes.bot <- genes.data[genes.data$Y == 0.5,]
  
  # Make Plotting Variables
  p.max <- -log10(p.min)
  y.max <-  max(round.up(p.max, decimals = 0), 8)
  if(!is.na(Gene) | !is.na(Region)[1]){
    x.min <- (gene.start - basepairs)
    x.max <- (gene.end + basepairs)
  } else{
    x.min <- (snp.pos - basepairs)
    x.max <- (snp.pos + basepairs)
  }
  
  y_offset <- y.max / 20
  
  if(!is.na(SecondarySNP)){
    x_offset = abs(x.max - x.min) / 150 * 15
    secondary.snp.pos <- new.results.data[new.results.data[,"snps"] == SecondarySNP, "pos"]
    secondary.snp.p <- -log10(new.results.data[new.results.data[,"snps"] == SecondarySNP, "p"])
    if(secondary.snp.pos < snp.pos){
      secondary.label.x.offset <- secondary.snp.pos - x_offset
      secondary.line.x.offset <- secondary.snp.pos - (x_offset / 3)
    } else {
      secondary.label.x.offset <- secondary.snp.pos + x_offset
      secondary.line.x.offset <- secondary.snp.pos + (x_offset / 3)
    }
    secondary.data <- new.results.data[new.results.data$pos > (secondary.snp.pos - (abs(x.max - x.min) / 6)) & new.results.data$pos < (secondary.snp.pos + (abs(x.max - x.min) / 6)), ]
    secondary.label.y.offset <- -log10(min(secondary.data[, "p"])) * 1.03
    if(abs(Nominal - secondary.label.y.offset) <= 1){
      secondary.label.y.offset <- max(secondary.label.y.offset, Nominal) * 1.03
    }
  }
  
  # Make Plot
  jpeg(width = 150, height = 160, units = "mm", res = 300, file = File.Name)
  layout(matrix(c(1, 2, 3), byrow = TRUE), heights = c(2, 6, 4))
  
  ## plot 1 - where SNPs are
  par(mar = c(0.6, 4, 4, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = rep(1, times = length(new.results.data.plot[,"snps"])), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted\nSNPs", las = 2, main = Plot.Title, xlim = c(x.min, x.max))
  
  ## plot 2 - actual manhattan p-value
  par(mar = c(0.5, 4, 0.6, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = -log10(new.results.data.plot[,"p"]), ylim = c(0, y.max), pch = 20, col = as.character(new.results.data.plot[,"Colour"]), xlab = "", ylab = expression(-log[10](italic(P))), cex = 1.5, xaxt = "n", xlim = c(x.min, x.max))
  points(x = snp.pos, y = snp.p, pch = 18, cex = 2, col = "#7D26CD")
  abline(h = Nominal, col = "blue", lty = "dashed")
  abline(h = Significant, col = "red", lty = "dashed")
  text(x = (snp.pos), y = (p.max + y_offset), labels = SNP)
  
  if(!is.na(SecondarySNP)){
    text(x = secondary.label.x.offset, y = secondary.label.y.offset, labels = SecondarySNP, cex = 0.8)
    segments(x0 = secondary.snp.pos, x1 = secondary.line.x.offset, y0 = secondary.snp.p, y1 = secondary.label.y.offset)
  }
  
  legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2"), col = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), fill = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), border = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), pt.cex = 2, cex = 1.2, bg = "white", box.lwd = 0, title = expression("r"^2))
  
  ## plot 3 - where the genes are
  par(mar = c(4, 4, 0.5, 4), mgp = c(2, 1, 0))
  plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", Chr), ylab="", xlim = c(x.min, x.max), ylim = c(0,3))
  # plot first line of genes
  odd = 1
  for(gene in 1:length(genes.top[,"Gene"])){
    lines(x = c(genes.top[gene,"Start"], genes.top[gene,"End"]), y = c(genes.top[gene,"Y"], genes.top[gene,"Y"]), lwd = 3, col = "#000080")
    if(length(genes.data[,"Gene"]) >= 10){
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
  # plot middle line of genes
  odd = 1
  for(gene in 1:length(genes.mid[,"Gene"])){
    lines(x = c(genes.mid[gene,"Start"], genes.mid[gene,"End"]), y = c(genes.mid[gene,"Y"], genes.mid[gene,"Y"]), lwd = 3, col = "#000080")
    if(length(genes.data[,"Gene"]) >= 10){
      length = abs(genes.mid[gene,"Start"] - genes.mid[gene,"End"])
      if(length > 8000 & odd%%2 == 0){ ## using the %% tells you if a number is odd/even odd == 1, even == 0
        text(x = (genes.mid[gene,"Start"] + genes.mid[gene,"End"])/2, y = genes.mid[gene, "Y"] - 0.2, labels = genes.mid[gene, "Gene"], cex = 0.8)
      } else{
        if(length > 8000 & odd%%2 == 1){
          text(x = (genes.mid[gene,"Start"] + genes.mid[gene,"End"])/2, y = genes.mid[gene, "Y"] + 0.2, labels = genes.mid[gene, "Gene"], cex = 0.8)
        }
      }
    } else{
      text(x = (genes.mid[gene,"Start"] + genes.mid[gene,"End"])/2, y = genes.mid[gene, "Y"] + 0.2, labels = genes.mid[gene, "Gene"])
    }
    odd = odd + 1
  }
  # plot third line of genes
  odd = 1
  for(gene in 1:length(genes.bot[,"Gene"])){
    lines(x = c(genes.bot[gene,"Start"], genes.bot[gene,"End"]), y = c(genes.bot[gene,"Y"], genes.bot[gene,"Y"]), lwd = 3, col = "#000080")
    if(length(genes.data[,"Gene"]) >= 10){
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

rm(Genes.Data, gene, genes.data, LD.colours, ld.data, LD.File, new.ld.data, new.results.data, new.results.data.plot, results.data, BP, Chr, CHR, basepairs, Nominal, P, p.max, p.min, Plot.Title, Significant, SNP, SNP.List, snp.p, snp.pos, x.max, x.min, y, y.max, y_offset)
  
  if(!is.na(SecondarySNP)){
    rm(x_offset, secondary.snp.pos, secondary.snp.p, secondary.label.x.offset, secondary.line.x.offset, secondary.data, secondary.label.y.offset)
  }

}
