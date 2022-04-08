#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Switch TOPBOT alleles
//
// \code{topbotChange} changes genotype calls for SNPs corresponding
// to the BOT naming convention used by Illumina. SNPs with genotypes
// of 0 will now have genotypes of 2 and vis-versa. Het calls remain
// unchanged.
//
// @param @param Genotypes An integer matrix of genotypes where each column represents a sample
// and each row represents a SNP. Genotypes are coded as -1 for missing, 0 for homozygous
// reference, 1 for heterozygous and 2 for homozygous alternative.
// @param topbot A character vector of Ts and Bs corresponding to TOPs and BOTs
// [[Rcpp::export]]
IntegerMatrix  topbotChange(IntegerMatrix Genotypes, CharacterVector topbot){
  const int noSNPs = Genotypes.nrow();
  const int noInds = Genotypes.ncol();
  IntegerMatrix newGenotypes(noSNPs,noInds);
  newGenotypes = Genotypes;
  for(int i = 0; i < noSNPs; ++i){
    if(topbot[i] == "B"){
      for(int j = 0; j < noInds; ++j){
        if(Genotypes(i,j) == 0 ) newGenotypes(i,j) = 3;
        if(Genotypes(i,j) == 2 ) newGenotypes(i,j) = 0;
        if(Genotypes(i,j) == 3 ) newGenotypes(i,j) = 2;
        if(Genotypes(i,j) == 1 ) newGenotypes(i,j) = 1;
      }
    }
  }
  return newGenotypes;
}
