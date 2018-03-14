# LocusZooms
## Make LocusZoom-like plots with your own LD matrix.

This script creates an R function to create regional Manhattan plots with points coloured according to LD and genes annotated beneath. Three example input files are included for test purposes, along with an example .jpg output.

  - Example.assoc.linear: A file of PLINK association results (only the SNP, BP, and P columns are essential)
  - Example.ld: A file of the LD between the SNP to be labelled (top-hit / SNP of interest) and the SNPs included in the PLINK results file
    - this file MUST have a column called "SNP_B" (containing a list of all the SNPs in the results file) and a column called "R2" (containing the R^2 LD value of each SNP). The SNP names MUST match the names in the SNP column of the results file.
 - Example.genes: A file of the genes within the region for use in the annotation step. This file must have 3 columns, Gene, Start, End. The UCSC_GRCh37_Genes_UniqueList.txt file can be used as this file.

### Example locus.zoom run:

```
# load necessary files into R
Example.assoc.linear <- read.delim("Example.assoc.linear", stringsAsFactors = FALSE, header = TRUE)
Example.ld <- read.table("Example.ld", stringsAsFactors = FALSE, header = TRUE)
UCSC_GRCh37_Genes_UniqueList.txt <- read.delim("UCSC_GRCh37_Genes_UniqueList.txt", stringsAsFactors = FALSE, header = TRUE)

# load the locuszoom function into R
source("functions/locus_zoom.R")

# create a LocusZoom-like plot
locus.zoom(CHR = Example.assoc.linear$CHR,                                 # the chromosome column in the results file 
           BP = Example.assoc.linear$BP,                                   # the position column in the results file
           P = Example.assoc.linear$P,                                     # the p-value column in the results file
           SNP.List = Example.assoc.linear$SNP,                            # the SNP names column in the resuls file
           SNP = "rs1008400",                                              # the SNP to be labelled in the plot
           LD.File = Example.ld,                                           # a file with LD values relevant to the SNP specified above
           Genes.Data = UCSC_GRCh37_Genes_UniqueList.txt,                  # a file of all the genes in the region / genome
           Plot.Title = "Association of FTO with BMI in Europeans",        # the plot title
           File.Name = "Example.jpg")                                      # the name of the file to save the plot to
```

One of SNP, Gene, or Region must be specified to create the plot:

 - SNP: specify the SNP to be annotated
 - Gene: specify the Gene to make the plot around (you will need to find out which SNP is going to be labelled and make an LD file relevant to that)
 - Region: specify the chromsome region you want to plot (must be specified as `c(chr, start, end)`


Other oppitional conditions are also available:

 - basepairs: specify how far either side of the SNP, Gene, or Region you want the plot to extend (defaults to 200000)
 - NonCoding: when using the UCSC gene list you can specify whether you want to plot the non-coding genes (defaults to FALSE)
 - Nominal: specify the nominal significance level to draw on the plot (in -log[10](_P_), default is 6 or _P_ = 1e-6)
 - Significant: specify the significance level to draw on the plot (in -log[10](_P_), default is 7.3 or _P_ = 5e-8) 
 - SecondarySNP: provide the ID of a second SNP to be labelled (will also include an arrow)
 
 
