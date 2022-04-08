#!/bin/bash

CHROM=20
threads=16


dir=Input/chr$CHROM
mkdir -p $dir

phased_data=../step1_compute_relatedness/1000GP.chr$CHROM.MAF5.vcf.gz



# 0. Unrelated set of samples
unrelated_samples=../step2_groups/Grouping/100unrelated.txt
out_unrelated=$dir/unrelated.chr$CHROM.vcf.gz
bcftools view --threads $threads -S $unrelated_samples -Oz -o $out_unrelated $phased_data
bcftools index --threads $threads -c $out_unrelated




# 1. Related set of samples
related_samples=../step2_groups/Grouping/Related_samples.txt
out_related=$dir/related.chr$CHROM.vcf.gz

bcftools view --threads $threads --force-samples -S $related_samples -Oz -o $out_related $phased_data
bcftools index --threads $threads -c $out_related





# 2. Parent-of-Origin inference
mkdir -p Out/chr$CHROM
groups=../step2_groups/Grouping/Full_callset.txt
maps=./POO_original/maps/chr$CHROM.b37.gmap.gz

./POO_original/bin/pooer -I $out_related -H $out_unrelated -M $maps -G $groups -R $CHROM -T $threads -O Out/chr$CHROM/chr$CHROM.prob.gz



