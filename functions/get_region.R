# Function to make regions out of variant or gene information
get.region <- function(snp.dat, snp, gene.dat, gene) {
  # If SNP is given:
  if (!is.na(snp)) {
    snp.ind = which(snp.dat$SNP == snp)
    
    if(length(snp.ind) == 0){
      stop(paste0("Your SNP (", snp, ") is not present in your data."))
    }
    
    snp.chr = snp.dat$CHR[snp.ind]
    snp.pos = snp.dat$BP[snp.ind]
    region = c(snp.chr, snp.pos, snp.pos)
  } else if (!is.na(gene)) {
    # If Gene is given
    gene.ind = which(gene.dat$Gene == gene)
    
    if(length(gene.ind) == 0){
      stop(paste0("Your gene (", gene, ") is not present in your data."))
    }
    
    gene.chr = gene.dat$Chrom[gene.ind]
    gene.start = gene.dat$Start[gene.ind]
    gene.end = gene.dat$End[gene.ind]
    region = c(gene.chr, gene.start, gene.end)
  } else {
    # If nothing was given
    stop("You must specify a SNP, Gene, or Region to plot")
  }
  return(region)
}