#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Calculate Missingness Proportions
//
// \code{calculateMissingness} is a function used to calculate the proportion of
// missing genotype calls for either a SNP or a sample
//
// @param genotypes An integer matrix of genotypes where each column represents a sample
// and each row represents a SNP. Genotypes are coded as -1 for missing, 0 for homozygous
// reference, 1 for heterozygous and 2 for homozygous alternative.
// [[Rcpp::export]]
NumericVector calculateMissingness(IntegerMatrix genotypes) {
  NumericVector proportion_missing(genotypes.ncol());
  int number_snps = genotypes.nrow();
  int number_isolates = genotypes.ncol();
  double number_snps_1 = genotypes.nrow();

  for (int i = 0; i < number_isolates; i++) {
    double number_missing = 0.0;
    for (int j = 0; j < number_snps; j++) {
      if(genotypes(j,i) == -1 ) number_missing += 1;
    }
    proportion_missing[i] = number_missing/number_snps_1;
  }
  return proportion_missing;
}
