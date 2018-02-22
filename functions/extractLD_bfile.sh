#!/bin/sh

# Usage: bash extractLD_bfile.sh plink_path bfile_prefix keep_file chr pos_start pos_end top_snp out_file
plink_path=$1
bfile_prefix=$2
keep_file=$3
chr=$4
pos_start=$5
pos_end=$6
top_snp=$7
out_file=$8

${plink_path}/plink \
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
    --out ld_files/${out_file}
