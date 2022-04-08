#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Group Combinations for Analysis
//
// Creates a data frame containing family IDs and individual IDs for each pair to be analysed. Each row
// corresponds to a unique pair.
//
// @param group A character vector of all family IDs
// [[Rcpp::export]]
CharacterMatrix groupPairs(CharacterVector group) {
  int number_isolates = group.size();
  CharacterMatrix group_pairs(number_isolates*(number_isolates - 1)/2,2);
  int k = 0;

  for (int i = 0; i < (number_isolates - 1); i++) {
    for(int j = (i + 1); j < number_isolates; j++){
      group_pairs(k,0) = group[i];
      group_pairs(k,1) = group[j];
      k += 1;
    }
  }
  return group_pairs;
}
