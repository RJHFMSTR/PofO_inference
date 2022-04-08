#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// IBD Segment Label
//
// \code{IBDLabel} is a function used to label unique IBD segments in a pair of samples for a particular chromosome by determining
// breakpoints in IBD vs non-IBD regions. IBD segments are labelled in sequential order genome wide.
//
// @param snp_id A numeric vector of SNP identifiers for IBD segments on the chromosome of interest.
// @param number_snps integer. The number of IBD SNPs for the chromosome of interest.
// [[Rcpp::export]]
IntegerVector IBDLabel(IntegerVector snp_id, const int number_snps){
  IntegerVector ibd_number(number_snps);

  ibd_number[0] = 1;
  for(int t = 1; t < number_snps; ++t){
    if(snp_id[t] - snp_id[t-1] == 1) ibd_number[t] = ibd_number[t-1];
    if(snp_id[t] - snp_id[t-1] > 1)  ibd_number[t] = ibd_number[t-1] + 1;
  }
  return ibd_number;
}
