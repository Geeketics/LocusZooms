# Function to subset relevant region from a dataset
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