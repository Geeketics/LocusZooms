######################
# Locus Zoom, Own LD Data
# Sptember 2016
# Tanya Flynn
# Uni of Otago

### Important Running Notes:
## You need to give the LD.Matrix a header & column called "snps" before using - MAKE SURE YOU HAVE YOUR HEADER LABELS IN THE RIGHT ORDER!!
## Genes.Data expects a data.frame with headers Gene, Start, End - this will be used to annotate the graph.

## LD Colours
# 1.0 - 0.8 = #FF0000
# 0.8 - 0.6 = #FFA500
# 0.6 - 0.4 = #00FF00
# 0.4 - 0.2 = #87CEFA
# 0.2 - 0.0 = #000080
# NA = #7F7F7F
# top-hit = #7D26CD

locus.zoom <- function(BP = NULL, P = NULL, SNP.List = NULL, SNP = NA, LD.Matrix = NULL, Chr = NULL, Genes.Data = NULL, Plot.Title = NULL, Nominal = 6, Significant = 7.3, File.Name = NULL){
  # Load Data
  results.data <- data.frame(snps = SNP.List, pos = BP, p = P)
  results.data[,"snps"] <- as.character(results.data[,"snps"])
  ld.data <- LD.Matrix
  genes.data <- Genes.Data

  # Extract Relevant LD
  p.min <- min(results.data[,"p"], na.rm = TRUE)
  if(is.na(SNP)){
    SNP <- results.data[results.data[,"p"] == p.min & !is.na(results.data[,"p"]), "snps"]
  }
  new.ld.data <- ld.data[,c("snps", SNP)]

  # Add LD to Results
  round.up <- function(x, decimals = 1){
    round(x + (5 * 10 ^ (-decimals - 1)), digits = decimals)
  }

  new.results.data <- merge(results.data, new.ld.data, by = "snps", all.x = TRUE)
  new.results.data[,"plot.ld"] <- round.up(new.results.data[,SNP], decimals = 1)
  new.results.data[new.results.data[,"plot.ld"] > 1 & !is.na(new.results.data[,"plot.ld"]), "plot.ld"] = 1
  
  LD.colours <- data.frame(LD = c(seq(from = 0, to = 1, by = 0.1), "NaN"), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2), "#7F7F7F"))
  
  new.results.data.plot <- merge(new.results.data, LD.colours, by.x = "plot.ld", by.y = "LD", all.x = TRUE)
  
  # Make Plotting Variables
  p.max <- -log10(p.min)
  ymax <-  max(round.up(p.max, decimals = 0), 8)
  x.min <-  min(new.results.data.plot[,"pos"], na.rm = TRUE)
  x.max <-  max(new.results.data.plot[,"pos"], na.rm = TRUE)
  snp.pos <- new.results.data.plot[new.results.data.plot[,"snps"] == SNP, "pos"]
  snp.p <- -log10(results.data[results.data[,"snps"] == SNP, "p"])
  
  y = rep(c(1.5, 0.5), times = length(genes.data[,"Gene"]))
  genes.data[,"Y"] = y[1:length(genes.data[,"Gene"])]
  
  # Make Plot
  jpeg(width = 1800, height = 1800, res = 300, file = File.Name)
  layout(matrix(c(1, 2, 3), byrow = TRUE), heights = c(2, 6, 3))
  
  ## plot 1 - where SNPs are
  par(mar = c(0.6, 4, 4, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = rep(1, times = length(new.results.data.plot[,"snps"])), axes = FALSE, pch = "|", xlab = "", ylab = "Plotted SNPs", las = 2, main = Plot.Title)
  
  ## plot 2 - actual manhattan p-value
  par(mar = c(0.5, 4, 0.6, 4), mgp = c(2, 1, 0))
  plot(x = new.results.data.plot[,"pos"], y = -log10(new.results.data.plot[,"p"]), ylim = c(0, ymax), pch = 20, col = as.character(new.results.data.plot[,"Colour"]), xlab = "", ylab = expression(-log[10](italic(P))), cex = 1.5, xaxt = "n")
  points(x = snp.pos, y = snp.p, pch = 18, cex = 2, col = "#7D26CD")
  text(x = snp.pos, y = (p.max * 1.2), labels = SNP)
  abline(h = Nominal, col = "blue", lty = "dashed")
  abline(h = Significant, col = "red", lty = "dashed")
  legend(x = "topright", legend = c("1.0", "0.8", "0.6", "0.4", "0.2"), col = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), fill = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), border = c("#FF0000", "#FFA500", "#00FF00", "#87CEFA", "#000080"), pt.cex = 2, cex = 1.2, bg = "white", box.lwd = 0)
  
  ## plot 3 - where the genes are
  par(mar = c(4, 4, 0.5, 4), mgp = c(2, 1, 0))
  plot(1, type = "n", yaxt = "n", xlab = paste("Position on Chromosome", Chr), ylab="", xlim = c(x.min, x.max), ylim = c(0,2))
  for(gene in 1:length(genes.data[,"Gene"])){
    lines(x = c(genes.data[gene,"Start"], genes.data[gene,"End"]), y = c(genes.data[gene,"Y"], genes.data[gene,"Y"]), lwd = 3, col = "#000080")
    text(x = (genes.data[gene,"Start"] + genes.data[gene,"End"])/2, y = genes.data[gene, "Y"] + 0.2, labels = genes.data[gene, "Gene"])
  }
  dev.off()
}
