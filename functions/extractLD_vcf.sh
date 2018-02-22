#!/bin/sh

# Usage: bash extractLD_vcf.sh plink_path vcf chr pos_start pos_end top_snp out_file

plink_path=$1
vcf=$2
chr=$3
pos_start=$4
pos_end=$5
top_snp=$6
out_file=$7

${plink_path}/plink \
    --vcf ${vcf} \
    --chr ${chr} \
    --from-bp ${pos_start} \
    --to-bp ${pos_end} \
    --r2 \
    --ld-window 100000000 \
    --ld-window-kb 100000000 \
    --ld-window-r2 0 \
    --ld-snp ${top_snp} \
    --out ld_files/${out_file}
