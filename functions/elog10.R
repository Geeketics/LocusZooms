# Function to convert string P-value into logged P:
elog10 <- function(p) {
  if (is.character(p) & grepl('e-', p)) {
    split_p <- base::strsplit(p, split = "e-")
    tmp <- unlist(split_p)
    res <- as.numeric(tmp[2]) - log10(as.numeric(tmp[1]))
  } else {
    res <- -log10(as.numeric(p))
  }
  return(res)
}