#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Calculate population allele frequencies from a matrix of genotypes.
//
// @param genotypes an integer matrix of genotypes where each column represents a sample
// and each row represents a SNP. Genotypes are coded as -1 for missing, 0 for homozygous
// reference, 1 for heterozygous and 2 for homozygous alternative.
// @param chromosomes a vector of chromosomes in the same order as the SNPs apear
// in \code{genotypes}. One chromosome per SNP.
// @param genders a vector of sample genders (1 = male, 2 = female)
//  in the same order as the samples apear in \code{genotypes}.
// [[Rcpp::export]]
NumericVector calculatePopAlleleFreq(IntegerMatrix genotypes, IntegerVector chromosomes, IntegerVector genders) {
  NumericVector pop_allele_freqs(genotypes.nrow());
  int number_samples = genotypes.ncol();
  int number_snps = genotypes.nrow();

  for (int t = 0; t < number_snps; t++) {
    if (chromosomes[t] != 23) {
      double A = 0, B = 0;
      for (int i = 0; i < number_samples; i++) {
        if (genotypes(t,i) == 0)  A += 2;
        if (genotypes(t,i) == 1){ A += 1; B += 1; }
        if (genotypes(t,i) == 2)  B += 2;
      }
      if(A + B == 0) pop_allele_freqs[t] = -1; else pop_allele_freqs[t] = A/(A+B);
    }
    if (chromosomes[t] == 23) {
      double A = 0, B = 0;
      for (int i = 0; i < number_samples; i++) {
        if (genotypes(t,i) == 0 && genders[i] == 1)  A += 1;
        if (genotypes(t,i) == 0 && genders[i] == 2)  A += 2;
        if (genotypes(t,i) == 1 && genders[i] == 2){ A += 1; B += 1; }
        if (genotypes(t,i) == 2 && genders[i] == 1)  B += 1;
        if (genotypes(t,i) == 2 && genders[i] == 2)  B += 2;
      }
      if (A + B == 0) pop_allele_freqs[t] = -1; else pop_allele_freqs[t] = A/(A+B);
    }
  }

  return pop_allele_freqs;
}
