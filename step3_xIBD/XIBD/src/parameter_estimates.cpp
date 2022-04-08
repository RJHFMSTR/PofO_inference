#include <Rcpp.h>
using namespace Rcpp;

//@ Internal Function
//
// Parameter Estimates
//
// Below are the functions used to perform the method-of-moments parameter
// estimation as in PLINK, Purcell et al. (2007). This method solves
// Ax=b for x in order to estimate omega values. The matrices A and b are
// created in this script.


//----- 2 haploid chromosomes, b and A -----

// The vector b in the equation Ax=b for 2 haploid chromosomes
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorHH(IntegerMatrix genotypes){
  IntegerVector Z(2);
  int IBS_0, IBS_1 = 0;
  int number_snps = genotypes.nrow();
  int number_snps_1 = 0;

  for(int i = 0; i < number_snps; ++i){
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      if(genotypes(i,0) == genotypes(i,1)){
        IBS_1 += 1;
      }
      number_snps_1 += 1;
    }
  }
  IBS_0 = number_snps_1 - IBS_1;
  Z[0] = IBS_0;
  Z[1] = IBS_1;
  return(Z);
}

// The matrix A in the equation Ax=b for 2 haploid chromosomes
// @param pop_allele_freqs A numeric vector of population allele frequencies for each SNP
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
NumericMatrix AmatrixHH(NumericVector pop_allele_freqs,IntegerMatrix genotypes){
NumericMatrix A(2,2);
int number_snps = pop_allele_freqs.size();
double i0_z0 = 0, i0_z1 = 0, i1_z0 = 0, i1_z1 = 0;
double p, q;

for(int i = 0; i < number_snps; ++i){
if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
p = pop_allele_freqs[i];
q = 1 - p;
i0_z0 += 2*p*q;
i1_z0 += pow(p,2)+pow(q,2);
i1_z1 += 1;
}
}
A(0,0) = i0_z0;
A(0,1) = i1_z0;
A(1,0) = i0_z1;
A(1,1) = i1_z1;

return(A);
}


//----- 1 haploid chromosome and 1 diploid chromosome, b and A -----

// The vector b in the equation Ax=b for 1 haploid chromosome and 1 diploid chromosome
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorHD(IntegerMatrix genotypes){
IntegerVector Z(2);
int IBS_0, IBS_1 = 0;
int number_snps = genotypes.nrow();
int number_snps_1 = 0;

for(int i = 0; i < number_snps; ++i){
if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
if(genotypes(i,0) == genotypes(i,1) || genotypes(i,0) == 1 || genotypes(i,1) == 1){
IBS_1 += 1;
}
number_snps_1 += 1;
}
}
IBS_0 = number_snps_1 - IBS_1;
Z[0] = IBS_0;
Z[1] = IBS_1;
return(Z);
}

// The matrix A in the equation Ax=b for 1 haploid and 1 diploid chromosome
// @param pop_allele_freqs A numeric vector of population allele frequencies for each SNP
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
NumericMatrix AmatrixHD(NumericVector pop_allele_freqs,IntegerMatrix genotypes){
  NumericMatrix A(2,2);
  int number_snps = pop_allele_freqs.size();
  double i0_z0 = 0, i0_z1 = 0, i1_z0 = 0, i1_z1 = 0;
  double p, q;

  for(int i = 0; i < number_snps; ++i){
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      p = pop_allele_freqs[i];
      q = 1 - p;
      i0_z0 += pow(p,2)*q + p*pow(q,2);
      i1_z0 += pow(p,3) + pow(q,3) + 2*pow(p,2)*q + 2*p*pow(q,2);
      i1_z1 += 1;
    }
  }
  A(0,0) = i0_z0;
  A(0,1) = i1_z0;
  A(1,0) = i0_z1;
  A(1,1) = i1_z1;

  return(A);
}


//----- 2 diploid chromosomes, b and A -----

// The vector b in the equation Ax=b for 2 diploid chromosomes
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorDD(IntegerMatrix genotypes){
  IntegerVector Z(3);
  int IBS_0 = 0, IBS_1 = 0, IBS_2 = 0;
  int number_snps = genotypes.nrow();

  for(int i = 0; i < number_snps; ++i){
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      if(genotypes(i,0) == genotypes(i,1)){
        IBS_2 += 1;
      }
      if((genotypes(i,0) == 1 && genotypes(i,1) == 0) || (genotypes(i,0) == 1 && genotypes(i,1) == 2) ||
         (genotypes(i,0) == 0 && genotypes(i,1) == 1) || (genotypes(i,0) == 2 && genotypes(i,1) == 1)){
        IBS_1 += 1;
      }
      if((genotypes(i,0) == 0 && genotypes(i,1) == 2) || (genotypes(i,0) == 2 && genotypes(i,1) == 0)){
        IBS_0 += 1;
      }
    }
  }
  Z[0] = IBS_0;
  Z[1] = IBS_1;
  Z[2] = IBS_2;
  return(Z);
}

// The matrix A in the equation Ax=b for 2 diploid chromosomes
// @param pop_allele_freqs A numeric vector of population allele frequencies for each SNP
// @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
NumericMatrix AmatrixDD(NumericVector pop_allele_freqs,IntegerMatrix genotypes){
  NumericMatrix A(3,3);
  int number_snps = pop_allele_freqs.size();
  double i0_z0 = 0, i0_z1 = 0, i0_z2 = 0, i1_z0 = 0, i1_z1 = 0, i1_z2 = 0;
  double i2_z0 = 0, i2_z1 = 0, i2_z2 = 0;
  double p, q;

  for(int i = 0; i < number_snps; ++i){
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      p = pop_allele_freqs[i];
      q = 1 - p;
      i0_z0 += 2*pow(p,2)*pow(q,2);
      i1_z0 += 4*pow(p,3)*q + 4*p*pow(q,3);
      i2_z0 += pow(p,4) + pow(q,4) + 4*pow(p,2)*pow(q,2);
      i1_z1 += 2*pow(p,2)*q + 2*p*pow(q,2);
      i2_z1 += pow(p,3) + pow(q,3) + pow(p,2)*q + p*pow(q,2);
      i2_z2 += 1;
    }
  }
  A(0,0) = i0_z0;
  A(0,1) = i1_z0;
  A(0,2) = i2_z0;
  A(1,0) = i0_z1;
  A(1,1) = i1_z1;
  A(1,2) = i2_z1;
  A(2,0) = i0_z2;
  A(2,1) = i1_z2;
  A(2,2) = i2_z2;
  return(A);
}

