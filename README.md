# LocusZooms
## Make LocusZoom-like plots with your own LD matrix.

This script creates an R function to create regional Manhattan plots with points coloured according to LD and genes annotated beneath. Three example input files are included for test purposes, along with an example .jpg output.

  - Example.assoc.linear: A file of PLINK association results (only the SNP, BP, and P columns are essential)
  - Example.ld: A file of the LD between the SNP to be labelled (top-hit / SNP of interest) and the SNPs included in the PLINK results file
    - this file MUST have a column called "SNP_B" (containing a list of all the SNPs in the results file) and a column called "R2" (containing the R^2 LD value of each SNP). The SNP names MUST match the names in the SNP column of the results file.
    - this file can be created for you by the locus_zoom.R script IF you have access to the Biochem servers and have rsIDs in your results file
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
locus.zoom(data = Example.assoc.linear,                                    # a data.frame (or a list of data.frames) with the columns CHR, BP, SNP, and P
           snp = "rs1008400",                                              # the SNP to be labelled in the plot
           ld.file = Example.ld,                                           # a file with LD values relevant to the SNP specified above
           genes.data = UCSC_GRCh37_Genes_UniqueList.txt,                  # a file of all the genes in the region / genome
           plot.title = "Association of FTO with BMI in Europeans",        # the plot title
           file.name = "Example.jpg")                                      # the name of the file to save the plot to
```

One of `snp`, `gene`, or `region` must be specified to create the plot:

 - snp: specify the SNP to be annotated
 - gene: specify the Gene to make the plot around
 - region: specify the chromsome region you want to plot (must be specified as `c(chr, start, end)`


Other oppitional conditions are also available:

 - `offset_bp`: specify how far either side of the `snp`, `gene`, or `region` you want the plot to extend (defaults to 200000)
 - `non-coding`: when using the UCSC gene list you can specify whether you want to plot the non-coding genes (defaults to FALSE)
 - `nominal`: specify the nominal significance level to draw on the plot (in -log[10](_P_), default is 6 or _P_ = 1e-6)
 - `significant`: specify the significance level to draw on the plot (in -log[10](_P_), default is 7.3 or _P_ = 5e-8) 
 - `secondary.snp`: provide the list of secondary SNP IDs (must match IDs in results file) to be highlighted on the plot
 - `secondary.label`: specify whether to label the secondary SNPs on the plot (defaults to FALSE)
 - `population`: specify the 1000 genomes population to use when calculating LD if ld.file = NULL (defaults to "EUR", options are "AFR", "AMR", "EAS", "EUR", and "SAS")
 - `sig.type`: specify whether the y-axis should be labelled as -log10(P) or -log10(BF) (defaults to "P", options are "P" or "BF")
 - `nplots`: specify whether multiple results plots will be saved into your jpeg file (e.g. plot two GWAS results one above another; defaults to FALSE)
 - `ignore.lead`: specify whether to ignore the SNP with the smallest P and use the SNP specified by 'snp' to centre the plot (defaults to FALSE)
 - `rsid.check`: specify whether to check if the SNPs are labelled with rsIDs - should only matter if script is calculating LD for you (defaults to TRUE)
 
