#! /bin/bash

# Script to generate LD of all the lead variants from plink loci output

LOCI=$1
ANCESTRY=$2

grep -v "^$" ${LOCI} | tr -s ' ' '\t' | sed -e 's/^\t//g' -e 's/^23/X/g' | cut -f1,3,4 | awk 'NR > 1 {print $2, $1":"$3 - 500000"-"$3 + 500000}' | tr ' ' '\t' > ${LOCI}.query

# Generate LD for each lead variant:
run_ld() {
	line=$1
	ancestry=$2
	rsid=$(echo ${line} | cut -d ' ' -f1)
	query=$(echo ${line} | cut -d ' ' -f2)
	chr=$(echo ${query} | cut -d ':' -f1)
	bcftools view -r ${query} -S /Volumes/archive/merrimanlab/riku/1kgp_sample_files/gazal_sample_list/${ancestry}_sample_list.txt /Volumes/archive/merrimanlab/reference_files/VCF/1000Genomes_vcf_files/Phase3_March2017/ALL.chr${chr}.*.vcf.gz -Oz -o results/locus_zooms/ld/${ancestry}/${rsid}.vcf.gz
	plink1.9b4.9 --vcf results/locus_zooms/ld/${ancestry}/${rsid}.vcf.gz --allow-no-sex --snps-only --r2 inter-chr --ld-snp ${rsid} --ld-window-r2 0 --out results/locus_zooms/ld/${ancestry}/${rsid}
}

export -f run_ld
cat ${LOCI}.query | parallel run_ld {} ${ANCESTRY}
