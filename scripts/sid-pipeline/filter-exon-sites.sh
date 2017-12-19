echo "filtering non-exons snps..."
parallel -u -t "$@" "../ensembl_exons.py < snps_{}.csv > exon_snps_{}.csv" ::: {1..19} X Y MT
