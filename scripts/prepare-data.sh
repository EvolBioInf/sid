#!/bin/bash

chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT"

# split input BAM into chromosomes (here for mouse)
parallel -u -t -j3 "samtools view -b C57BL_6NJ.bam {} > C57BL_6NJ_{}.bam" ::: $chromosomes

# remove duplicate reads
for c in $chromosomes; do
    java -jar picard.jar MarkDuplicates INPUT=C57BL_6NJ_$c.bam OUTPUT=C57BL_6NJ_$c.dedup.bam METRICS_FILE=dedup-metrics-$c.txt REMOVE_DUPLICATES=true
done

# create pileups
parallel -u -t -j3 "samtools mpileup -C50 -q1 C57BL_6NJ_{}.dedup.bam | gzip -c > C57BL_6NJ_{}.plp.gz" ::: $chromosomes
