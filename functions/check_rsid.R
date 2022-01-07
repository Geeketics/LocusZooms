# Function to check if the input variants have rsIDs
check.rsid <- function(snp = NULL) {
  # Stop if CHR:POS ID:
  if (all(!grepl('rs', snp))) {
    stop("Your SNP column does not have rsIDs")
  }
  # Stop if there are duplicate SNPs:
  if (length(which(duplicated(snp))) > 0) {
    stop("There are duplicate rsIDs in your results file - Please remove them before running again")
  }
}