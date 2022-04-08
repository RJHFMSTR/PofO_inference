#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Get Sample Pairs
//
// \code{isolatePairs} is a function used to get all combinations of
// pairs from a given family ID vector and sample ID vector
//
// @param fid A character vector of family IDs
// @param iid A character vector of individual IDs
// [[Rcpp::export]]
CharacterMatrix isolatePairs(CharacterVector fid, CharacterVector iid) {
  int number_isolates = fid.size();
  CharacterMatrix isolate_pairs(number_isolates*(number_isolates - 1)/2,4);
  int k = 0;

  for (int i = 0; i < (number_isolates - 1); i++) {
    for(int j = (i + 1); j < number_isolates; j++){
      isolate_pairs(k,0) = fid[i];
      isolate_pairs(k,1) = iid[i];
      isolate_pairs(k,2) = fid[j];
      isolate_pairs(k,3) = iid[j];
      k += 1;
    }
  }
  return isolate_pairs;
}
