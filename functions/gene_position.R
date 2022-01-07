# Function to plot the gene tracks and labels properly:
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