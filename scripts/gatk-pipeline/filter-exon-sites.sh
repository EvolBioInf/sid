echo "filtering non-exons snps..."
parallel -u -t "$@" "../vcfsite2csv.sh < het_snps_{}.vcf | ../ensembl_exons.py > exon_snps_{}.csv" ::: {1..19} X Y MT
