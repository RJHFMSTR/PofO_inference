#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Get Genotypes from Haplotypes
//
// \code{haplotypeToGenotype} is a function used to extract genotypes from haplotypes, where a
// haplotype consists of two alleles and a genotype is a single representation of these alleles.
//
// @param haplotypes An integer matrix of haplotypes where each row represents a sample
// and two adjacent columns represent the alleles for a single SNP.
// @param chromosomes An integer vector of chromosomes; one chromosome per SNP.
// @param genders An integer vector of genders (1 or 2) of the samples in the haplotypes matrix.
// [[Rcpp::export]]
IntegerMatrix haplotypeToGenotype(IntegerMatrix haplotypes, IntegerVector chromosomes, IntegerVector genders){
  IntegerMatrix genotypes(haplotypes.ncol()/2, haplotypes.nrow());
  int number_samples = haplotypes.nrow();
  int number_snps = haplotypes.ncol();

  for (int i = 0; i < number_samples; i++) {
    for (int j = 0; j < number_snps; j=j+2) {
      if(haplotypes(i,j) == 1 && haplotypes(i,j+1) == 1) genotypes(j/2,i) = 0;
      if(haplotypes(i,j) == 2 && haplotypes(i,j+1) == 2) genotypes(j/2,i) = 2;
      if(haplotypes(i,j) == 0 && haplotypes(i,j+1) == 0) genotypes(j/2,i) = -1;
      if(chromosomes[j] == 23 & genders[i] == 1 & ((haplotypes(i,j) == 1 && haplotypes(i,j+1) == 2) || (haplotypes(i,j) == 2 && haplotypes(i,j+1) == 1))) genotypes(j/2,i) = -1;
      if(chromosomes[j] == 23 & genders[i] == 2 & ((haplotypes(i,j) == 1 && haplotypes(i,j+1) == 2) || (haplotypes(i,j) == 2 && haplotypes(i,j+1) == 1))) genotypes(j/2,i) = 1;
      if(chromosomes[j] != 23 & ((haplotypes(i,j) == 1 && haplotypes(i,j+1) == 2) || (haplotypes(i,j) == 2 && haplotypes(i,j+1) == 1))) genotypes(j/2,i) = 1;
    }
  }
  return genotypes;
}
