#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// Calculate Haplotype Frequencies
//
// This script contains 3 functions used in the calculation of
// haplotype frequencies between pairs of SNPs. SNPs to condition
// on are already assumed to have been found.


//----- round decimal -----

// Round digits to specified decimal places
// @param number A number to round
// @param number The number of digits to round to
double round_decimal_h(double number, int digits){
  double new_number;
  new_number = ceil( number * pow(10,digits) - 0.4999999 )/ pow(10,digits);
  if(new_number == -0) return 0;
  else return new_number;
}


//----- population allele frequency -----

// Calculate population allele frequency for a single SNP
// @param genotypes An integer matrix of genotypes with number of rows equal to number of SNPs
// and number of columns equal to the number of samples. Genotypes are coded as -1 for missing,
// 0 for homozygous reference, 1 for heterozygous and 2 for homozygous alternative.
// @param chromosome An interger value from 1 - 23 giving the chromosome of the SNP
// @param genders An integer vector of the genders for samples in genotypes matrix
double calculatePopAlleleFreqSingle(IntegerVector genotypes, const int chromosome, IntegerVector genders) {
  double pop_allele_freqs;
  int number_samples = genotypes.size();

  if (chromosome != 23) {
    double A = 0, B = 0;
    for (int i = 0; i < number_samples; i++) {
      if (genotypes[i] == 0)  A += 2;
      if (genotypes[i] == 1){ A += 1; B += 1; }
      if (genotypes[i] == 2)  B += 2;
    }
    if(A + B == 0) pop_allele_freqs = -1; else pop_allele_freqs = A/(A+B);
  }
  if (chromosome == 23) {
    double A = 0, B = 0;
    for (int i = 0; i < number_samples; i++) {
      if (genotypes[i] == 0 && genders[i] == 1)  A += 1;
      if (genotypes[i] == 0 && genders[i] == 2)  A += 2;
      if (genotypes[i] == 1 && genders[i] == 2){ A += 1; B += 1; }
      if (genotypes[i] == 2 && genders[i] == 1)  B += 1;
      if (genotypes[i] == 2 && genders[i] == 2)  B += 2;
    }
    if (A + B == 0) pop_allele_freqs = -1; else pop_allele_freqs = A/(A+B);
  }

  return pop_allele_freqs;
}


//----- haplotype frequency for a single SNP -----

// Calculate haplotype frequency for a single SNP
// @param genotypes An integer matrix of genotypes with number of rows equal to number of SNPs
// and number of columns equal to the number of samples. Genotypes are coded as -1 for missing,
// 0 for homozygous reference, 1 for heterozygous and 2 for homozygous alternative.
// @param genders An integer vector of the genders for samples in genotypes matrix
// @param current_snp Integer denoting the numeric ID of the SNP of interest
// @param condition_snp Integer denoting the numeric ID of the SNP to condition on
// @param chromosome An interger value from 1 - 23 giving the chromosome of the SNP
// @param polyrootR An R function to calculate the roots of a polynomial
NumericMatrix haplotypeFreq(IntegerMatrix genotypes, IntegerVector genders, const int current_snp, const int condition_snp, const int chromosome, Function polyrootR ){
  const int number_samples = genotypes.ncol();
  double p11, p12, p21, p22 ;
  NumericMatrix condition_matrix(1,6) ;
  IntegerVector missing_geno(number_samples);

  int a1 = 0, b1 = 0, c1 = 0, d1 = 0, e1 = 0, f1 = 0, g1 = 0, h1 = 0, i1 = 0;
  int a2 = 0, b2 = 0, c2 = 0, d2 = 0;
  int X, Y, Z, Q, find_roots = 0;
  double coef0, coef1, coef2, coef3;
  double A1 = 0, B1 = 0, C1 = 0, D1 = 0, A, B, C, D;
  ComplexVector roots(3);

  for(int ind = 0; ind < number_samples; ++ind){
    if((chromosome != 23) || (chromosome == 23 && genders[ind] == 2)){
      find_roots = 1;
      if(genotypes(condition_snp,ind) == 2 && genotypes(current_snp,ind) == 2) a1 += 1;
      if(genotypes(condition_snp,ind) == 2 && genotypes(current_snp,ind) == 1) b1 += 1;
      if(genotypes(condition_snp,ind) == 2 && genotypes(current_snp,ind) == 0) c1 += 1;
      if(genotypes(condition_snp,ind) == 1 && genotypes(current_snp,ind) == 2) d1 += 1;
      if(genotypes(condition_snp,ind) == 1 && genotypes(current_snp,ind) == 1) e1 += 1;
      if(genotypes(condition_snp,ind) == 1 && genotypes(current_snp,ind) == 0) f1 += 1;
      if(genotypes(condition_snp,ind) == 0 && genotypes(current_snp,ind) == 2) g1 += 1;
      if(genotypes(condition_snp,ind) == 0 && genotypes(current_snp,ind) == 1) h1 += 1;
      if(genotypes(condition_snp,ind) == 0 && genotypes(current_snp,ind) == 0) i1 += 1;
    }
    if(chromosome == 23 && genders[ind] == 1){
      if(genotypes(condition_snp,ind) == 2 && genotypes(current_snp,ind) == 2) a2 += 1;
      if(genotypes(condition_snp,ind) == 2 && genotypes(current_snp,ind) == 0) b2 += 1;
      if(genotypes(condition_snp,ind) == 0 && genotypes(current_snp,ind) == 2) c2 += 1;
      if(genotypes(condition_snp,ind) == 0 && genotypes(current_snp,ind) == 0) d2 += 1;
    }
    if(genotypes(current_snp,ind) == -1) missing_geno[ind] = 1;
    if(genotypes(current_snp,ind) != -1) missing_geno[ind] = 0;
  }
  if(find_roots == 1){
    X = 2*a1 + b1 + d1;
    Y = 2*i1 + h1 + f1;
    Z = 2*c1 + b1 + f1;
    Q = 2*g1 + h1 + d1;

    coef0 = -X*Y;
    coef1 = Z*Q + e1*(Z+Q) - e1*(X+Y) + pow(e1,2) + X*Y;
    coef2 = -e1*(Z+Q) - 2*pow(e1,2) - pow(e1,2) + e1*(X+Y);
    coef3 = 2*pow(e1,2);

    NumericVector polycoef = NumericVector::create (coef0,coef1,coef2,coef3) ;
    roots = polyrootR(polycoef) ; // R function to get roots of poynomial. Imaginary part ignored.

    double p = 0.5, q = 0.5;
    double p0, q0;
    double like;
    for(int R = 0; R < 3; R++){
      p0 = roots[R].r;
      p0 = round_decimal_h(p0,4);
      q0 = 1 - p0;
      if(R == 0){
        like = p0/q0;
        p = p0;
        q = q0;
      }
      if(R == 1 | R == 2){
        if(p0/q0 > like){
          like = p0/q0;
          p = p0;
          q = q0;
        }
      }
    }
    A1 = 2*a1 + b1 + d1 + p*e1;
    B1 = 2*c1 + b1 + f1 + q*e1;
    C1 = 2*g1 + h1 + d1 + q*e1;
    D1 = 2*i1 + h1 + f1 + p*e1;
  }
  A = A1 + a2;
  B = B1 + b2;
  C = C1 + c2;
  D = D1 + d2;

  double sumABCD = A + B + C + D;

  p11 = A/sumABCD ;
  p12 = B/sumABCD ;
  p21 = C/sumABCD ;
  p22 = D/sumABCD ;

  // set missing genotypes at current SNP to missing at condition SNP
  IntegerVector genotypes_temp(number_samples) ;
  for(int ind = 0; ind < number_samples; ++ind){
    if(missing_geno[ind] == 1) genotypes_temp[ind] = -1;
    if(missing_geno[ind] != 1) genotypes_temp[ind] = genotypes(condition_snp,ind);
  }

  condition_matrix(0,0) = condition_snp + 1;
  condition_matrix(0,1) = p11; // pba;
  condition_matrix(0,2) = p12; // pbA;
  condition_matrix(0,3) = p21; // pBa;
  condition_matrix(0,4) = p22; // pBA;
  condition_matrix(0,5) = calculatePopAlleleFreqSingle(genotypes_temp, chromosome, genders);

  return condition_matrix;
}


//----- haplotype frequencies for all SNPs -----

// Calculate haplotype frequencues for all SNPs
// @param genotypes An integer matrix of genotypes with number of rows equal to number of SNPs
// and number of columns equal to the number of samples. Genotypes are coded as -1 for missing,
// 0 for homozygous reference, 1 for heterozygous and 2 for homozygous alternative.
// @param condition_snps_0 An integer matrix containing the numeric IDs of each SNP and numeric IDs
// of their condition SNP. The first SNP of each chromosome has condition_snp=-1 since no conditioning
// is performed on this SNP. Numeric IDs of the SNPs are seq(1,number.snps,1).
// @param genders An integer vector of the genders for samples in genotypes matrix
// @param chromosome An interger value from 1 - 23 giving the chromosome of the SNP
// @param polyrootR An R function to calculate the roots of a polynomial (polyroot)
// [[Rcpp::export]]
NumericMatrix calculateHaplotypeFreq(IntegerMatrix genotypes, IntegerMatrix condition_snps_0, IntegerVector genders, const int chromosome, Function polyrootR){
  int number_snps = genotypes.nrow();
  int current_snp;
  int condition_snp;
  NumericMatrix condition_snps(number_snps,6);
  NumericMatrix condition_info(1,6);

  // get condition info for first SNP
  condition_snps(0,0) = -1 ;
  condition_snps(0,1) = 0 ;
  condition_snps(0,2) = 0 ;
  condition_snps(0,3) = 0 ;
  condition_snps(0,4) = 0 ;
  condition_snps(0,5) = calculatePopAlleleFreqSingle(genotypes(0,_), chromosome, genders);

  // get condition info for remaining SNPs
  for(int snp = 1; snp < number_snps; ++snp){
    current_snp = condition_snps_0(snp,0) - 1 ;
    condition_snp = condition_snps_0(snp,1) - 1 ;
    condition_info = haplotypeFreq(genotypes, genders, current_snp, condition_snp, chromosome, polyrootR) ;
    condition_snps(snp,0) = condition_info(0,0) ;
    condition_snps(snp,1) = condition_info(0,1) ;
    condition_snps(snp,2) = condition_info(0,2) ;
    condition_snps(snp,3) = condition_info(0,3) ;
    condition_snps(snp,4) = condition_info(0,4) ;
    condition_snps(snp,5) = condition_info(0,5) ;
  }
  return condition_snps ;
}
