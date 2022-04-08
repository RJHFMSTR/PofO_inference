#!/bin/bash

#SBATCH --job-name Gr_groups
#SBATCH --array=1-1
#SBATCH --output Err/%x_%A-%a.out
#SBATCH --error Err/%x_%A-%a.err
#SBATCH --partition normal
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mail-type ALL
#SBATCH --export=NONE
#SBATCH --mem 40G
#SBATCH --time 1-0



CHROM=20

threads=16


VCF=1000GP.chr20.MAF5.vcf.gz
mkdir -p Plink
plink --vcf $VCF --threads 16 --make-bed --out Plink/chr$CHROM

mkdir -p King

/users/rhofmeis/commands/KING/king -b Plink/chr$CHROM.bed --related --degree 3 --cpus $threads --prefix King/chr$CHROM










