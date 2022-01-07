# Function to get the LD information of specified population from the 1000 Genomes data (March 2017 release):
# NOTE: the input SNP MUST be in rsID format, not CHR:POS-based.
# NOTE: This function will leave/save the LD information in the working directory for future reference (e.g. if the user wanted to use the same LD information)
get.ld <- function(region, snp, population) {
  ld.snp = snp
  
  if (region[1] == "23") {
    region[1] = "X"
  }
  
  vcf.filename = "POP_chrZZ.no_relatives.no_indel.biallelic.vcf.gz"
  vcf.filename = gsub(pattern = 'ZZ', replacement = region[1], vcf.filename)
  if(population == "TAMA"){
    vcf.filename = gsub(pattern = 'POP', replacement = "AFR_AMR_EAS_EUR", vcf.filename)
  } else{
    vcf.filename = gsub(pattern = 'POP', replacement = population, vcf.filename)
  }
  
  # check necessary 1000 genomes file can be reached
  if(!(vcf.filename %in% list.files(path = paste0("/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/", population)))){
    stop(paste0("The file ", paste0("/Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/", population, "/", vcf.filename), " cannot be found."))
  }
  
  # gsub the command and filename for chr, start/end positions and the population:
  base.command = "source ~/.bashrc;
  bcftools view \
    --regions ZZ:Y1-Y2 \
    --output-type z \
    --output-file tmp.vcf.gz \
    /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/POP/1000VCF;

  plink \
    --vcf tmp.vcf.gz \
    --allow-no-sex \
    --snps-only \
    --r2 \
    --inter-chr \
    --ld-snp SNP \
    --ld-window-r2 0 \
    --out POP_region_ZZ.Y1-Y2_SNP;

  rm tmp.vcf.gz POP_region_ZZ.Y1-Y2_SNP.nosex"
  
  base.command = gsub(pattern = "\n ", replacement = "", base.command)
  base.command = gsub(pattern = 'ZZ', replacement = region[1], base.command)
  base.command = gsub(pattern = 'Y1', replacement = region[2], base.command)
  base.command = gsub(pattern = 'Y2', replacement = region[3], base.command)
  base.command = gsub(pattern = "1000VCF", replacement = vcf.filename, base.command)
  base.command = gsub(pattern = 'POP', replacement = population, base.command)
  base.command = gsub(pattern = 'SNP', replacement = ld.snp, base.command)
  
  # Make a system call to run the bcftools/plink command.
  # I'm only assigning it to a variable to suppress any possible form of
  # messages/outputs from the command, just in case
  # ignores any errors when running the LD command, but will output the error to your screen
  messages = system(base.command, ignore.stdout = TRUE, intern = TRUE, ignore.stderr = TRUE)
  
  # Import the LD data:
  ld.file = "POP_region_ZZ.Y1-Y2_SNP.ld"
  ld.file = gsub(pattern = 'ZZ', replacement = region[1], ld.file)
  ld.file = gsub(pattern = 'Y1', replacement = region[2], ld.file)
  ld.file = gsub(pattern = 'Y2', replacement = region[3], ld.file)
  ld.file = gsub(pattern = 'POP', replacement = population, ld.file)
  ld.file = gsub(pattern = 'SNP', replacement = ld.snp, ld.file)
  
  # Check the ld file was made & import it
  if(ld.file %in% list.files(pattern = ".ld")){
    ld = read.table(ld.file, stringsAsFactors = FALSE, header = TRUE)
  } else{
    message("Top SNP / specified SNP not in 1000 Genomes biallelic SNPs")
    ld = data.frame(CHR_A = NA, BP_A = NA, SNP_A = NA, CHR_B = NA, BP_B = NA, SNP_B = NA, R2 = NA)
  }
  # return the ld file data
  return(ld)
}