echo "filtering synonymous mutations..."
parallel -u -t "$@" "../nonsynonymous.py exon_snps_{}.csv ../sid-local-output/raw_{}.csv.gz > nonsynonymous_exon_snps_{}.csv" ::: {1..19} X Y MT

