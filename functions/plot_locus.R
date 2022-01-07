# Function to create LocusZoom style plot (without gene track):
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