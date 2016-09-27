# LocusZooms
## Make LocusZoom-like plots with your own LD matrix.

This script creates an R function to create regional Manhattan plots with points coloured according to LD and genes annotated beneath. Three example input files are included for test purposes, along with an example .jpg output.

 - Example.assoc.linear: A file of PLINK association results (only the SNP, BP, and P columns are essential)
 - Example.ld: A matrix of the LD between the SNPs included in the PLINK results file - this file MUST have SNP names for the headers and a column called "snps" of the SNP names. The SNP names MUST match the names in the SNP column of the results file.
 - Example.genes: A file of the genes within the region for use in the annotation step. This file must have 3 columns, Gene, Start, End.

### Example locus.zoom run:
```
source(locus_zoom.R)
locus.zoom(BP = Example.assoc.linear$BP, P = Example.assoc.linear$P, SNP.List = Example.assoc.linear$SNP, LD.Matrix = Example.ld, Chr = 16, Genes.Data = Example.genes, Plot.Title = "Association of FTO with BMI in West Polynesians", File.Name = "Example.jpg")
```
Three oppitional conditions are also available:

 - SNP: specify the SNP to be annotated/LD to be plotted relevant to (default is the most significant SNP)
 - Nominal: specify the nominal significance level (default is 1e-6)
 - Significant: specify the significance level (default is 5e-8) 

