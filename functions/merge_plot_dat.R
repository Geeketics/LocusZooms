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