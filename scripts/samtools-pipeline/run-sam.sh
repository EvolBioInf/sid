#!/bin/bash

chrom=$1
bam="$HOME/data/C57BL_6NJ_$chrom.dedup.bam"
raw="raw_$chrom.vcf"
called="called_$chrom.vcf"
het="het_snps_$chrom.vcf"

ref=$HOME/data/Mus_musculus.GRCm38.68.dna_sm.toplevel.fa

# parameters from ftp-mouse.sanger.ac.uk/current_snps/README
samtools mpileup -t DP,DV,DP4,SP,DPR,INFO/DPR -E -Q 0 -pm3 -F0.25 -d500 -v -f $ref $bam > $raw &&
bcftools call -mv -f GQ,GP -p 0.99 $raw > $called &&
grep -v INDEL $called | grep -e "0/1" -e "1/2" -e "0/2" > $het
