#!/bin/sh

bfile_prefix=$1
keep_file=$2
chr=$3
pos_start=$4
pos_end=$5
top_snp=$6
file_prefix=$7

/Users/tanyaflynn/Executables/plink-1.90b3.29/plink \
    --bfile ${bfile_prefix} \
    --keep ${keep_file} \
    --chr ${chr} \
    --from-bp ${pos_start} \
    --to-bp ${pos_end} \
    --r2 \
    --ld-window 100000000 \
    --ld-window-kb 100000000 \
    --ld-window-r2 0 \
    --ld-snp ${top_snp} \
    --out ld_files/${file_prefix}
