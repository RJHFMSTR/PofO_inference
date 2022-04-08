#!/bin/bash
CHROM=20

threads=16
VCF=../step1_compute_relatedness/1000GP.chr20.MAF5.vcf.gz
scaffold=chr20.scaffold.vcf.gz

# Step 5.1  produce phasing tree

mkdir -p Bingraph/


maps=./shapeit4/maps/chr$CHROM.b37.gmap.gz
bingraph=Bingraph/chr$CHROM.graph


./shapeit4/bin/shapeit4.2 -I $VCF -R $CHROM -S $scaffold -M $maps -T $threads --bingraph $bingraph


# Step 5.2 sample 1,000 random haplotypes

mkdir -p Sampling

./shapeit4/tools/bingraphsample/bin/bingraphsample -I $bingraph --thread $threads --collapse 1000 -O Sampling/Bin.sampling.chr$CHROM.vcf.gz


# Step 5.3 identify segment of high confidence from the 1,000 phased haplotypes
python src/uncertainty_segments.py



# Step 5.4 output the most likely phased halotypes

./shapeit4/tools/bingraphsample/bin/bingraphsample -I $bingraph --solve -O Phasing_best_guess.vcf.gz





