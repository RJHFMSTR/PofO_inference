#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Binary Matrix of IBD
//
// Creates a binary matrix of IBD (1) and non-IBD (0) with each row representing a single SNP
// and each column representing a unique pair. The number of rows is equal to the total
// number of SNPs and the number of columns is equal to the number of pairs.
//
// @param chromosomes A character vector containing the corresponding chromosome for each SNP
// @param positions_bp A numeric vector containing the corresponding bp position for each SNP
// @param number_pairs Numeric. The total number of pairs analysed
// @param ibd_pairs_colnumbers A numeric vector corresponding to column numbers in the output matrix
// where each unique number refers to a unique pair with IBD inferred
// @param ibd_chromosomes A character vector containing the chromosome for each detected IBD segment
// @param ibd_start_bp A numeric vector containing the base-pair position for the start of each detected IBD segment
// @param ibd_stop_bp A numeric vector containing the base-pair position for the end of each detected IBD segment
// [[Rcpp::export]]
IntegerMatrix IBDMatrix(CharacterVector chromosomes, NumericVector positions_bp, int number_pairs, IntegerVector ibd_pairs_colnumbers,
            CharacterVector ibd_chromosomes, NumericVector ibd_start_bp, NumericVector ibd_stop_bp) {
  IntegerMatrix ibd_m(chromosomes.size(), number_pairs);

  for(int i=0; i<ibd_chromosomes.size(); i++){
    for(int j=0; j<chromosomes.size(); j++){
      if(chromosomes[j] == ibd_chromosomes[i] && positions_bp[j] >= ibd_start_bp[i] && positions_bp[j] <= ibd_stop_bp[i]){
        ibd_m(j,ibd_pairs_colnumbers[i]-1) = 1;
      }
    }
  }
  return ibd_m;
}
