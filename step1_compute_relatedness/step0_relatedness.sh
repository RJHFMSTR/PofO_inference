#!/bin/bash

CHROM=20

threads=16


VCF=1000GP.chr20.MAF5.vcf.gz
mkdir -p Plink
plink --vcf $VCF --threads 16 --make-bed --out Plink/chr$CHROM

mkdir -p King

./KING_software/king -b Plink/chr$CHROM.bed --related --degree 3 --cpus $threads --prefix King/chr$CHROM










