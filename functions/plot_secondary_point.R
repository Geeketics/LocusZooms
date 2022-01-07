# Function to plot secondary SNPs:
plot.secondary.point <- function(data, snps, lead.snp, plot.data, plot.var, label = FALSE, circle = TRUE) {
  
  lead.ind = which(data$SNP == lead.snp)
  lead.pos = data$BP[lead.ind]
  
  data = data[which(data$SNP != lead.snp), ]
  data = data[data$BP >= plot.var[2] & data$BP <= plot.var[3], ]
  
  
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