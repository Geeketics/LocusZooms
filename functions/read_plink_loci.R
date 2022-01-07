# Function to read in and pull out relevant info from PLINK clump output:
read.plink.loci <- function(file = NULL) {
  if (is.null(file)) {
    stop('You must provide a file for reading')
  }
  data = read.table(file, stringsAsFactors = FALSE, header = TRUE)
  data = data[,c(1,3:5)]
  return(data)
}