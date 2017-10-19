#!/bin/sh

vcf=$1
chr=$2
pos_start=$3
pos_end=$4
top_snp=$5
file_prefix=$6

/Users/tanyaflynn/Executables/plink-1.90b3.29/plink \
    --vcf ${vcf} \
    --chr ${chr} \
    --from-bp ${pos_start} \
    --to-bp ${pos_end} \
    --r2 \
    --ld-window 100000000 \
    --ld-window-kb 100000000 \
    --ld-window-r2 0 \
    --ld-snp ${top_snp} \
    --out ld_files/${file_prefix}
