#!/bin/bash

chrom=$1
input="$HOME/data/C57BL_6NJ_$chrom.dedup.bam"
raw="raw_output_$chrom.vcf"
snps="raw_snps_$chrom.vcf"
filtered="filtered_snps_$chrom.vcf"
het="het_snps_$chrom.vcf"

ref=$HOME/data/Mus_musculus.GRCm38.68.dna_sm.toplevel.fa

gatk -T HaplotypeCaller -R $ref -I $input -L $chrom --genotyping_mode DISCOVERY -stand_call_conf 20 -o $raw &&
gatk -T SelectVariants -R $ref -V $raw -selectType SNP -o $snps &&
gatk -T VariantFiltration -R $ref -V $snps --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filterName "snpfilter" -o $filtered &&
grep PASS $filtered | grep -e "0/1" -e "1/2" > $het

