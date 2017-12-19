#!/bin/bash
chrom=$1

sid_params="-R -m likelihood_ratio"
# sid_params="-R -m local"

# Using the following may speed up computation, if /dev/shm is a RAM-mounted filesystem (tmpfs)
# (should be the case when using Ubuntu)
# tmpdir=$(mktemp -d -p /dev/shm/sid)

inputgz="$HOME/data/C57BL_6NJ_${chrom}.plp.gz"
raw="raw_$chrom.csv.gz"
snps="snps_$chrom.csv"

zcat $inputgz > ${tmpdir:-.}/input.plp &&
    sid $sid_params ${tmpdir:-.}/input.plp | gzip -c > $raw &&
    zgrep ,het, $raw > $snps

rm -r ${tmpdir:-.}/input.plp

