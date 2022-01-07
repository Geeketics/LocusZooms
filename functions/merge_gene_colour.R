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