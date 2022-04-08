#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// HMM Model 2
//
// All the functions needed from the HMM, including emission probabilities
// transition probabilities, alpha, beta, gamma and viterbi. Haplotype
// frequencies are calculated in a different script.
// Model is based on Albretchsen et al. 2009


//----- round decimal -----

// Round digits to specified decimal places
// @param number A number to round
// @param number The number of digits to round to
double roundDecimal(double number, int digits){
  double new_number;
  new_number = ceil( number * pow(10,digits) - 0.4999999 )/ pow(10,digits);
  if(new_number == -0) return 0;
  else return new_number;
}


//----- emission probabilities -----

// The emission probabilities for 2 haploid chromosomes
// @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
// @param genotype_1 The genotype for isolate 1 from the pair for SNP i
// @param genotype_2 The genotype for isolate 2 from the pair for SNP i
// @param ibd The IBD state
double emissionProbHH(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  double alt_allele_freq = 1 - pop_allele_freq;
  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,2), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal( pop_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal( alt_allele_freq, 6);

  return 0;
}


// The emission probabilities for 1 haploid chromosome and 1 diploid chromosome
// @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
// @param genotype_1 The genotype for isolate 1 from the pair for SNP i
// @param genotype_2 The genotype for isolate 2 from the pair for SNP i
// @param ibd The IBD state
// @param male_column The haploid isolate from the pair. Either 1 or 2
// @param female_column The diploid isolate from the pair. Either 1 or 2
double emissionProbHD(double pop_allele_freq, int genotype_1, int genotype_2, int ibd, int male_column, int female_column) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  int geno_male, geno_female;

  if(male_column == 1) { geno_male = genotype_1; geno_female = genotype_2; }
  if(male_column == 2) { geno_male = genotype_2; geno_female = genotype_1; }

  if(geno_female == 0 && geno_male == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,3), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 0 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 0) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,3), 6);

  if(geno_female == 0 && geno_male == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 1) return 0;
  if(geno_female == 1 && geno_male == 0 && ibd == 1) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 1) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 1) return 0;
  if(geno_female == 2 && geno_male == 2 && ibd == 1) return roundDecimal( pow(alt_allele_freq,2), 6);

  return 0;
}


// The emission probabilities for 2 diploid chromosomes
// @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
// @param genotype_1 The genotype for isolate 1 from the pair for SNP i
// @param genotype_2 The genotype for isolate 2 from the pair for SNP i
// @param ibd The IBD state
double emissionProbDD(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,4), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 0) return roundDecimal( 4*pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,4), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,3), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 1) return roundDecimal( (pow(pop_allele_freq,2)*alt_allele_freq + pop_allele_freq*pow(alt_allele_freq,2)), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 1) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 1) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal( pow(alt_allele_freq,3), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 2) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 2) return roundDecimal( 2*pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 2) return roundDecimal( pow(alt_allele_freq,2), 6);

  return 0;
}


//----- transition probabilities -----

// The transition probabilities for 2 haploid chromosomes
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbHH(double omega_0, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}


// The transition probabilities for 1 haploid and 1 diploid chromosome
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbHD(double omega_0, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}


// The transition probabilities for 2 diploid chromosomes
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbDD(double omega_0, double omega_1, double omega_2, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  double T02 = (exp(-alpha * omega_1 * dist_M) * omega_2)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_0 * omega_1)/(omega_1 - 1) + omega_2;
  double T20 = (exp(-alpha * omega_1 * dist_M) * omega_0)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_2 * omega_1)/(omega_1 - 1) + omega_0;

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T02), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 0 && ibd_current == 2) return roundDecimal( T02, 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_0), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1 + exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 2) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_2), 6);
  if(ibd_previous == 2 && ibd_current == 0) return roundDecimal( T20, 6);
  if(ibd_previous == 2 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 2 && ibd_current == 2) return roundDecimal( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T20), 6);

  return 0;
}


//----- genotype error probabilities -----

// The genotyping error probability for 1 haploid chromosome
// @param truth The true genotype
// @param observed The observed genotype
// @param error The genotype error rate
double genotypeErrorH(int truth, int observed, double error){
  if(truth == observed) return 1 - error;
  else return error;
}


// The genotyping error probability for 1 diploid chromosome
// @param truth The true genotype
// @param observed The observed genotype
// @param error The genotype error rate
double genotypeErrorD(int truth, int observed, double error){
  if(truth == 0 && observed ==0) return pow(1 - error,2);
  if(truth == 0 && observed ==1) return 2 * (1 - error) * error;
  if(truth == 0 && observed ==2) return pow(error,2);
  if(truth == 1 && observed ==0) return (1 - error) * error;
  if(truth == 1 && observed ==1) return pow(1 - error,2) + pow(error,2);
  if(truth == 1 && observed ==2) return (1 - error) * error;
  if(truth == 2 && observed ==0) return pow(error,2);
  if(truth == 2 && observed ==1) return 2 * (1 - error) * error;
  if(truth == 2 && observed ==2) return pow(1 - error,2);

  return 0;
}


//----- joint probabilities -----

// Joint probabilities for 2 haploid chromosomes
// @param haplotypeFreq Numeric vecotor of haplotype frequencies for a SNP
// @param ibd The ibd state of SNP j
NumericVector jointProbHH(NumericVector haplotypeFreq, const int ibd){
  double pBA, pBa, pbA, pba ;
  NumericVector joint(16) ;
  const int chrom = 23 ;

  pBA = haplotypeFreq[3] ;
  pBa = haplotypeFreq[2] ;
  pbA = haplotypeFreq[1] ;
  pba = haplotypeFreq[0] ;

  if(ibd == 0) joint = NumericVector::create (pow(pBA,2), pBA*pBa, pBA*pBa, pow(pBa,2), pBA*pbA, pBA*pba, pBa*pbA, pBa*pba, pbA*pBA, pbA*pBa, pba*pBA, pba*pBa, pow(pbA,2), pbA*pba, pba*pbA, pow(pba,2)) ;
  if(ibd == 1) joint = NumericVector::create (pBA, 0, 0, pBa, 0, 0, 0, 0, 0, 0, 0, 0, pbA, 0, 0, pba) ;
  return joint ;
}


// Joint probabilities for 1 haploid and 1 diploid chromosome
// @param haplotypeFreq Numeric vecotor of haplotype frequencies for a SNP
// @param ibd The ibd state of SNP j
NumericVector jointProbHD(NumericVector haplotypeFreq, const int ibd){
  double pBA, pBa, pbA, pba ;
  NumericVector joint(36) ;
  const int chrom = 23 ;

  pBA = haplotypeFreq[3] ;
  pBa = haplotypeFreq[2] ;
  pbA = haplotypeFreq[1] ;
  pba = haplotypeFreq[0] ;

  if(ibd == 0){
    joint[0]=pow(pBA,3);          joint[1]=pow(pBA,2)*pBa;    joint[2]=2*pow(pBA,2)*pBa;                    joint[3]=2*pBA*pow(pBa,2);                    joint[4]=pBA*pow(pBa,2);    joint[5]=pow(pBa,3);
    joint[6]=pow(pBA,2)*pbA;      joint[7]=pow(pBA,2)*pba;    joint[8]=2*pBA*pbA*pBa;                       joint[9]=2*pBA*pBa*pba;                       joint[10]=pow(pBa,2)*pbA;   joint[11]=pow(pBa,2)*pba;
    joint[12]=2*pow(pBA,2)*pbA;   joint[13]=2*pBA*pbA*pBa;    joint[14]=2*pow(pBA,2)*pba + 2*pBa*pbA*pBA;   joint[15]=2*pbA*pow(pBa,2) + 2*pBA*pba*pBa;   joint[16]=2*pBa*pba*pBA;    joint[17]=2*pow(pBa,2)*pba;
    joint[18]=2*pBA*pow(pbA,2);   joint[19]=2*pBA*pbA*pba;    joint[20]=2*pBA*pba*pbA + 2*pBa*pow(pbA,2);   joint[21]=2*pBA*pow(pba,2) + 2*pBa*pbA*pba;   joint[22]=2*pBa*pba*pbA;    joint[23]=2*pBa*pow(pba,2);
    joint[24]=pow(pbA,2)*pBA;     joint[25]=pow(pbA,2)*pBa;   joint[26]=2*pbA*pba*pBA;                      joint[27]=2*pbA*pba*pBa;                      joint[28]=pow(pba,2)*pBA;   joint[29]=pow(pba,2)*pBa;
    joint[30]=pow(pbA,3);         joint[31]=pow(pbA,2)*pba;   joint[32]=2*pow(pbA,2)*pba;                   joint[33]=2*pow(pba,2)*pbA;                   joint[34]=pow(pba,2)*pbA;   joint[35]=pow(pba,3);
  }
  if(ibd == 1){
    joint[0]=pow(pBA,2);  joint[1]=0;   joint[2]=pBA*pBa;   joint[3]=pBA*pBa;   joint[4]=0;   joint[5]=pow(pBa,2);
    joint[6]=0;           joint[7]=0;   joint[8]=0;         joint[9]=0;         joint[10]=0;  joint[11]=0;
    joint[12]=pBA*pbA;    joint[13]=0;  joint[14]=pBA*pba;  joint[15]=pbA*pBa;  joint[16]=0;  joint[17]=pBa*pba;
    joint[18]=pBA*pbA;    joint[19]=0;  joint[20]=pBa*pbA;  joint[21]=pBA*pba;  joint[22]=0;  joint[23]=pBa*pba;
    joint[24]=0;          joint[25]=0;  joint[26]=0;        joint[27]=0;        joint[28]=0;  joint[29]=0;
    joint[30]=pow(pbA,2); joint[31]=0;  joint[32]=pbA*pba;  joint[33]=pba*pbA;  joint[34]=0;  joint[35]=pow(pba,2);
  }
  return joint ;
}


// Joint probabilities for 2 diploid chromosomes
// @param haplotypeFreq Numeric vecotor of haplotype frequencies for a SNP
// @param ibd The ibd state of SNP j
NumericVector jointProbDD(NumericVector haplotypeFreq, const int ibd){
  double pBA, pBa, pbA, pba ;
  NumericVector joint(81) ;
  const int chrom = 1 ;

  pBA = haplotypeFreq[3] ;
  pBa = haplotypeFreq[2] ;
  pbA = haplotypeFreq[1] ;
  pba = haplotypeFreq[0] ;


  if(ibd == 0){
    joint[0]=pow(pBA,4);                joint[1]=2*pow(pBA,3)*pBa;                          joint[2]=pow(pBA,2)*pow(pBa,2);   joint[3]=2*pow(pBA,3)*pBa;                          joint[4]=4*pow(pBA,2)*pow(pBa,2);                       joint[5]=2*pow(pBa,3)*pBA;                          joint[6]=pow(pBA,2)*pow(pBa,2);   joint[7]=2*pow(pBa,3)*pBA;                          joint[8]=pow(pBa,4);
    joint[9]=2*pow(pBA,3)*pbA;          joint[10]=2*pow(pBA,2)*pbA*pBa + 2*pow(pBA,3)*pba;  joint[11]=2*pow(pBA,2)*pba*pBa;   joint[12]=4*pow(pBA,2)*pbA*pBa;                     joint[13]=4*pow(pBA,2)*pba*pBa + 4*pBA*pow(pBa,2)*pbA;  joint[14]=4*pow(pBa,2)*pBA*pba;                     joint[15]=2*pow(pBa,2)*pBA*pbA;   joint[16]=2*pow(pBa,3)*pbA + 2*pow(pBa,2)*pBA*pba;  joint[17]=2*pow(pBa,3)*pba;
    joint[18]=pow(pBA,2)*pow(pbA,2);    joint[19]=2*pow(pBA,2)*pbA*pba;                     joint[20]=pow(pBA,2)*pow(pba,2);  joint[21]=2*pBA*pBa*pow(pbA,2);                     joint[22]=4*pba*pbA*pBA*pBa;                            joint[23]=2*pBA*pBa*pow(pba,2);                     joint[24]=pow(pBa,2)*pow(pbA,2);  joint[25]=2*pow(pBa,2)*pbA*pba;                     joint[26]=pow(pBa,2)*pow(pba,2);
    joint[27]=2*pow(pBA,3)*pbA;         joint[28]=4*pow(pBA,2)*pbA*pBa;                     joint[29]=2*pow(pBa,2)*pBA*pbA;   joint[30]=2*pow(pBA,3)*pba + 2*pow(pBA,2)*pbA*pBa;  joint[31]=4*pow(pBA,2)*pba*pBa + 4*pBA*pow(pBa,2)*pbA;  joint[32]=2*pow(pBa,2)*pBA*pba + 2*pow(pBa,3)*pbA;  joint[33]=2*pow(pBA,2)*pBa*pba;   joint[34]=4*pow(pBa,2)*pba*pBA;                     joint[35]=2*pow(pBa,3)*pba;
    joint[36]=4*pow(pBA,2)*pow(pbA,2);  joint[37]=4*pow(pBA,2)*pbA*pba + 4*pBA*pow(pbA,2)*pBa; joint[38]=4*pba*pbA*pBA*pBa;   joint[39]=4*pow(pBA,2)*pbA*pba + 4*pBA*pow(pbA,2)*pBa; joint[40]=4*pow(pBA,2)*pow(pba,2) + 8*pba*pbA*pBA*pBa + 4*pow(pBa,2)*pow(pbA,2); joint[41]=4*pow(pBa,2)*pba*pbA + 4*pBa*pow(pba,2)*pBA;    joint[42]=4*pba*pbA*pBA*pBa; joint[43]=4*pow(pBa,2)*pba*pbA + 4*pBa*pow(pba,2)*pBA; joint[44]=4*pow(pBa,2)*pow(pba,2);
    joint[45]=2*pBA*pow(pbA,3);         joint[46]=4*pow(pbA,2)*pBA*pba;                     joint[47]=2*pBA*pbA*pow(pba,2);   joint[48]=2*pow(pbA,2)*pBA*pba + 2*pow(pbA,3)*pBa;  joint[49]=4*pBA*pbA*pow(pba,2) + 4*pBa*pow(pbA,2)*pba;  joint[50]=2*pow(pba,2)*pBa*pbA + 2*pow(pba,3)*pBA;  joint[51]=2*pBa*pba*pow(pbA,2);   joint[52]=4*pBa*pbA*pow(pba,2);                     joint[53]=2*pow(pba,3)*pBa;
    joint[54]=pow(pBA,2)*pow(pbA,2);    joint[55]=2*pow(pbA,2)*pBA*pBa;                     joint[56]=pow(pbA,2)*pow(pBa,2);  joint[57]=2*pbA*pba*pow(pBA,2);                     joint[58]=4*pba*pbA*pBA*pBa;                            joint[59]=2*pbA*pba*pow(pBa,2);                     joint[60]=pow(pba,2)*pow(pBA,2);  joint[61]=2*pow(pba,2)*pBA*pBa;                     joint[62]=pow(pBa,2)*pow(pba,2);
    joint[63]=2*pBA*pow(pbA,3);         joint[64]=2*pow(pbA,2)*pBA*pba + 2*pow(pbA,3)*pBa;  joint[65]=2*pBa*pba*pow(pbA,2);   joint[66]=4*pow(pbA,2)*pBA*pba;                     joint[67]=4*pBA*pbA*pow(pba,2) + 4*pBa*pow(pbA,2)*pba;  joint[68]=4*pBa*pbA*pow(pba,2);                     joint[69]=2*pBA*pbA*pow(pba,2);   joint[70]=2*pow(pba,2)*pBa*pbA + 2*pow(pba,3)*pBA;  joint[71]=2*pow(pba,3)*pBa;
    joint[72]=pow(pbA,4);               joint[73]=2*pow(pbA,3)*pba;                         joint[74]=pow(pbA,2)*pow(pba,2);  joint[75]=2*pow(pbA,3)*pba;                         joint[76]=4*pow(pbA,2)*pow(pba,2);                      joint[77]=2*pow(pba,3)*pbA;                         joint[78]=pow(pbA,2)*pow(pba,2);  joint[79]=2*pow(pba,3)*pbA;                         joint[80]=pow(pba,4);
  }
  if(ibd == 1){
    joint[0]=pow(pBA,3);                        joint[1]=pow(pBA,2)*pBa;              joint[2]=0;   joint[3]=pow(pBA,2)*pBa;              joint[4]=pow(pBA,2)*pBa + pow(pBa,2)*pBA;   joint[5]=pow(pBa,2)*pBA;    joint[6]=0;   joint[7]=pow(pBa,2)*pBA;    joint[8]=pow(pBa,3);
    joint[9]=pow(pBA,2)*pbA;                    joint[10]=pow(pBA,2)*pba;             joint[11]=0;  joint[12]=pBA*pBa*pbA;                joint[13]=pBA*pba*pBa + pBA*pBa*pbA;        joint[14]=pBa*pBA*pba;      joint[15]=0;  joint[16]=pow(pBa,2)*pbA;   joint[17]=pow(pBa,2)*pba;
    joint[18]=0;                                joint[19]=0;                          joint[20]=0;  joint[21]=0;                          joint[22]=0;                                joint[23]=0;                joint[24]=0;  joint[25]=0;                joint[26]=0;
    joint[27]=pow(pBA,2)*pbA;                   joint[28]=pBA*pBa*pbA;                joint[29]=0;  joint[30]=pow(pBA,2)*pba;             joint[31]=pBA*pba*pBa + pBA*pBa*pbA;        joint[32]=pow(pBa,2)*pbA;   joint[33]=0;  joint[34]=pBa*pBA*pba;      joint[35]=pow(pBa,2)*pba;
    joint[36]=pBA*pow(pbA,2) + pow(pBA,2)*pbA;  joint[37]=pBA*pbA*pba + pBA*pbA*pBa;  joint[38]=0;  joint[39]=pBA*pbA*pba + pBA*pbA*pBa;  joint[40]=pow(pBA,2)*pba + pBA*pow(pba,2) + pow(pBa,2)*pbA + pBa*pow(pbA,2);          joint[41]=pBa*pba*pbA + pBa*pba*pBA;  joint[42]=0; joint[43]=pBa*pba*pbA + pBa*pba*pBA; joint[44]=pBa*pow(pba,2) + pow(pBa,2)*pba;
    joint[45]=pow(pbA,2)*pBA;                   joint[46]=pbA*pba*pBA;                joint[47]=0;  joint[48]=pow(pbA,2)*pBa;             joint[49]=pbA*pBa*pba + pbA*pba*pBA;        joint[50]=pow(pba,2)*pBA;   joint[51]=0;  joint[52]=pba*pbA*pBa;      joint[53]=pow(pba,2)*pBa;
    joint[54]=0;                                joint[55]=0;                          joint[56]=0;  joint[57]=0;                          joint[58]=0;                                joint[59]=0;                joint[60]=0;  joint[61]=0;                joint[62]=0;
    joint[63]=pow(pbA,2)*pBA;                   joint[64]=pow(pbA,2)*pBa;             joint[65]=0;  joint[66]=pbA*pba*pBA;                joint[67]=pbA*pBa*pba + pbA*pba*pBA;        joint[68]=pba*pbA*pBa,      joint[69]=0;  joint[70]=pow(pba,2)*pBA;   joint[71]=pow(pba,2)*pBa;
    joint[72]=pow(pbA,3);                       joint[73]=pow(pbA,2)*pba;             joint[74]=0;  joint[75]=pow(pbA,2)*pba;             joint[76]=pow(pbA,2)*pba + pow(pba,2)*pbA;  joint[77]=pow(pba,2)*pbA;   joint[78]=0;  joint[79]=pow(pba,2)*pbA;   joint[80]=pow(pba,3);
  }
  if(ibd == 2){
    joint[0]=pow(pBA,2);  joint[1]=0;   joint[2]=0;   joint[3]=0;   joint[4]=2*pBA*pBa;               joint[5]=0;   joint[6]=0;   joint[7]=0;   joint[8]=pow(pBa,2);
    joint[9]=0;           joint[10]=0;  joint[11]=0;  joint[12]=0;  joint[13]=0;                      joint[14]=0;  joint[15]=0;  joint[16]=0;  joint[17]=0;
    joint[18]=0;          joint[19]=0;  joint[20]=0;  joint[21]=0;  joint[22]=0;                      joint[23]=0;  joint[24]=0;  joint[25]=0;  joint[26]=0;
    joint[27]=0;          joint[28]=0;  joint[29]=0;  joint[30]=0;  joint[31]=0;                      joint[32]=0;  joint[33]=0;  joint[34]=0;  joint[35]=0;
    joint[36]=2*pBA*pbA;  joint[37]=0;  joint[38]=0;  joint[39]=0;  joint[40]=2*pBA*pba + 2*pBa*pbA;  joint[41]=0;  joint[42]=0;  joint[43]=0;  joint[44]=2*pBa*pba;
    joint[45]=0;          joint[46]=0;  joint[47]=0;  joint[48]=0;  joint[49]=0;                      joint[50]=0;  joint[51]=0;  joint[52]=0;  joint[53]=0;
    joint[54]=0;          joint[55]=0;  joint[56]=0;  joint[57]=0;  joint[58]=0;                      joint[59]=0;  joint[60]=0;  joint[61]=0;  joint[62]=0;
    joint[63]=0;          joint[64]=0;  joint[65]=0;  joint[66]=0;  joint[67]=0;                      joint[68]=0;  joint[69]=0;  joint[70]=0;  joint[71]=0;
    joint[72]=pow(pbA,2); joint[73]=0;  joint[74]=0;  joint[75]=0;  joint[76]=2*pbA*pba;              joint[77]=0;  joint[78]=0;  joint[79]=0;  joint[80]=pow(pba,2);
  }
  return joint ;
}


//----- true genotypes -----

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerVector GijC(int gender1, int gender2, int chrom){
  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    IntegerVector Gij(81) ;
    Gij[0]=0;   Gij[1]=0;   Gij[2]=0;   Gij[3]=0;   Gij[4]=0;   Gij[5]=0;   Gij[6]=0;   Gij[7]=0;   Gij[8]=0;   Gij[9]=0;   Gij[10]=0;  Gij[11]=0;  Gij[12]=0;  Gij[13]=0;
    Gij[14]=0;  Gij[15]=0;  Gij[16]=0;  Gij[17]=0;  Gij[18]=0;  Gij[19]=0;  Gij[20]=0;  Gij[21]=0;  Gij[22]=0;  Gij[23]=0;  Gij[24]=0;  Gij[25]=0;  Gij[26]=0;
    Gij[27]=1;  Gij[28]=1;  Gij[29]=1;  Gij[30]=1;  Gij[31]=1;  Gij[32]=1;  Gij[33]=1;  Gij[34]=1;  Gij[35]=1;  Gij[36]=1;  Gij[37]=1;  Gij[38]=1;  Gij[39]=1;  Gij[40]=1;
    Gij[41]=1;  Gij[42]=1;  Gij[43]=1;  Gij[44]=1;  Gij[45]=1;  Gij[46]=1;  Gij[47]=1;  Gij[48]=1;  Gij[49]=1;  Gij[50]=1;  Gij[51]=1;  Gij[52]=1;  Gij[53]=1;
    Gij[54]=2;  Gij[55]=2;  Gij[56]=2;  Gij[57]=2;  Gij[58]=2;  Gij[59]=2;  Gij[60]=2;  Gij[61]=2;  Gij[62]=2;  Gij[63]=2;  Gij[64]=2;  Gij[65]=2;  Gij[66]=2;  Gij[67]=2;
    Gij[68]=2;  Gij[69]=2;  Gij[70]=2;  Gij[71]=2;  Gij[72]=2;  Gij[73]=2;  Gij[74]=2;  Gij[75]=2;  Gij[76]=2;  Gij[77]=2;  Gij[78]=2;  Gij[79]=2;  Gij[80]=2;
    return Gij ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerVector Gij(16) ;
    Gij = IntegerVector::create (0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2) ;
    return Gij ;
  }
  if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
    IntegerVector Gij(36) ;
    Gij[0]=0;   Gij[1]=0;   Gij[2]=0;   Gij[3]=0;   Gij[4]=0;   Gij[5]=0;   Gij[6]=0;   Gij[7]=0;   Gij[8]=0;   Gij[9]=0;   Gij[10]=0;  Gij[11]=0;
    Gij[12]=1;  Gij[13]=1;  Gij[14]=1;  Gij[15]=1;  Gij[16]=1;  Gij[17]=1;  Gij[18]=1;  Gij[19]=1;  Gij[20]=1;  Gij[21]=1;  Gij[22]=1;  Gij[23]=1;
    Gij[24]=2;  Gij[25]=2;  Gij[26]=2;  Gij[27]=2;  Gij[28]=2;  Gij[29]=2;  Gij[30]=2;  Gij[31]=2;  Gij[32]=2;  Gij[33]=2;  Gij[34]=2;  Gij[35]=2;
    return Gij ;
  }
  return 0;
}

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerVector GikC(int gender1, int gender2, int chrom){
  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    IntegerVector Gik(81) ;
    Gik[0]=0;   Gik[1]=0;   Gik[2]=0;   Gik[3]=0;   Gik[4]=0;   Gik[5]=0;   Gik[6]=0;   Gik[7]=0;   Gik[8]=0;   Gik[9]=1;   Gik[10]=1;  Gik[11]=1;  Gik[12]=1;  Gik[13]=1;
    Gik[14]=1;  Gik[15]=1;  Gik[16]=1;  Gik[17]=1;  Gik[18]=2;  Gik[19]=2;  Gik[20]=2;  Gik[21]=2;  Gik[22]=2;  Gik[23]=2;  Gik[24]=2;  Gik[25]=2;  Gik[26]=2;
    Gik[27]=0;  Gik[28]=0;  Gik[29]=0;  Gik[30]=0;  Gik[31]=0;  Gik[32]=0;  Gik[33]=0;  Gik[34]=0;  Gik[35]=0;  Gik[36]=1;  Gik[37]=1;  Gik[38]=1;  Gik[39]=1;  Gik[40]=1;
    Gik[41]=1;  Gik[42]=1;  Gik[43]=1;  Gik[44]=1;  Gik[45]=2;  Gik[46]=2;  Gik[47]=2;  Gik[48]=2;  Gik[49]=2;  Gik[50]=2;  Gik[51]=2;  Gik[52]=2;  Gik[53]=2;
    Gik[54]=0;  Gik[55]=0;  Gik[56]=0;  Gik[57]=0;  Gik[58]=0;  Gik[59]=0;  Gik[60]=0;  Gik[61]=0;  Gik[62]=0;  Gik[63]=1;  Gik[64]=1;  Gik[65]=1;  Gik[66]=1;  Gik[67]=1;
    Gik[68]=1;  Gik[69]=1;  Gik[70]=1;  Gik[71]=1;  Gik[72]=2;  Gik[73]=2;  Gik[74]=2;  Gik[75]=2;  Gik[76]=2;  Gik[77]=2;  Gik[78]=2;  Gik[79]=2;  Gik[80]=2;
    return Gik ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerVector Gik(16) ;
    Gik = IntegerVector::create (0,0,0,0,2,2,2,2,0,0,0,0,2,2,2,2) ;
    return Gik ;
  }
  if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
    IntegerVector Gik(36) ;
    Gik[0]=0;   Gik[1]=0;   Gik[2]=0;   Gik[3]=0;   Gik[4]=0;   Gik[5]=0;   Gik[6]=2;   Gik[7]=2;   Gik[8]=2;   Gik[9]=2;   Gik[10]=2;  Gik[11]=2;
    Gik[12]=0;  Gik[13]=0;  Gik[14]=0;  Gik[15]=0;  Gik[16]=0;  Gik[17]=0;  Gik[18]=2;  Gik[19]=2;  Gik[20]=2;  Gik[21]=2;  Gik[22]=2;  Gik[23]=2;
    Gik[24]=0;  Gik[25]=0;  Gik[26]=0;  Gik[27]=0;  Gik[28]=0;  Gik[29]=0;  Gik[30]=2;  Gik[31]=2;  Gik[32]=2;  Gik[33]=2;  Gik[34]=2;  Gik[35]=2;
    return Gik ;
  }
  return 0;
}

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerVector GhjC(int gender1, int gender2, int chrom){
  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    IntegerVector Ghj(81) ;
    Ghj[0]=0;   Ghj[1]=0;   Ghj[2]=0;   Ghj[3]=1;   Ghj[4]=1;   Ghj[5]=1;   Ghj[6]=2;   Ghj[7]=2;   Ghj[8]=2;   Ghj[9]=0;   Ghj[10]=0;  Ghj[11]=0;  Ghj[12]=1;  Ghj[13]=1;
    Ghj[14]=1;  Ghj[15]=2;  Ghj[16]=2;  Ghj[17]=2;  Ghj[18]=0;  Ghj[19]=0;  Ghj[20]=0;  Ghj[21]=1;  Ghj[22]=1;  Ghj[23]=1;  Ghj[24]=2;  Ghj[25]=2;  Ghj[26]=2;
    Ghj[27]=0;  Ghj[28]=0;  Ghj[29]=0;  Ghj[30]=1;  Ghj[31]=1;  Ghj[32]=1;  Ghj[33]=2;  Ghj[34]=2;  Ghj[35]=2;  Ghj[36]=0;  Ghj[37]=0;  Ghj[38]=0;  Ghj[39]=1;  Ghj[40]=1;
    Ghj[41]=1;  Ghj[42]=2;  Ghj[43]=2;  Ghj[44]=2;  Ghj[45]=0;  Ghj[46]=0;  Ghj[47]=0;  Ghj[48]=1;  Ghj[49]=1;  Ghj[50]=1;  Ghj[51]=2;  Ghj[52]=2;  Ghj[53]=2;
    Ghj[54]=0;  Ghj[55]=0;  Ghj[56]=0;  Ghj[57]=1;  Ghj[58]=1;  Ghj[59]=1;  Ghj[60]=2;  Ghj[61]=2;  Ghj[62]=2;  Ghj[63]=0;  Ghj[64]=0;  Ghj[65]=0;  Ghj[66]=1;  Ghj[67]=1;
    Ghj[68]=1;  Ghj[69]=2;  Ghj[70]=2;  Ghj[71]=2;  Ghj[72]=0;  Ghj[73]=0;  Ghj[74]=0;  Ghj[75]=1;  Ghj[76]=1;  Ghj[77]=1;  Ghj[78]=2;  Ghj[79]=2;  Ghj[80]=2;
    return Ghj ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerVector Ghj(16) ;
    Ghj = IntegerVector::create (0,0,2,2,0,0,2,2,0,0,2,2,0,0,2,2) ;
    return Ghj ;
  }
  if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
    IntegerVector Ghj(36) ;
    Ghj[0]=0;   Ghj[1]=0;   Ghj[2]=1;   Ghj[3]=1;   Ghj[4]=2;   Ghj[5]=2;   Ghj[6]=0;   Ghj[7]=0;   Ghj[8]=1;   Ghj[9]=1;   Ghj[10]=2;  Ghj[11]=2;
    Ghj[12]=0;  Ghj[13]=0;  Ghj[14]=1;  Ghj[15]=1;  Ghj[16]=2;  Ghj[17]=2;  Ghj[18]=0;  Ghj[19]=0;  Ghj[20]=1;  Ghj[21]=1;  Ghj[22]=2;  Ghj[23]=2;
    Ghj[24]=0;  Ghj[25]=0;  Ghj[26]=1;  Ghj[27]=1;  Ghj[28]=2;  Ghj[29]=2;  Ghj[30]=0;  Ghj[31]=0;  Ghj[32]=1;  Ghj[33]=1;  Ghj[34]=2;  Ghj[35]=2;
    return Ghj ;
  }
  return 0;
}

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerVector GhkC(int gender1, int gender2, int chrom){
  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    IntegerVector Ghk(81) ;
    Ghk[0]=0;   Ghk[1]=1;   Ghk[2]=2;   Ghk[3]=0;   Ghk[4]=1;   Ghk[5]=2;   Ghk[6]=0;   Ghk[7]=1;   Ghk[8]=2;   Ghk[9]=0;   Ghk[10]=1;  Ghk[11]=2;  Ghk[12]=0;  Ghk[13]=1;
    Ghk[14]=2;  Ghk[15]=0;  Ghk[16]=1;  Ghk[17]=2;  Ghk[18]=0;  Ghk[19]=1;  Ghk[20]=2;  Ghk[21]=0;  Ghk[22]=1;  Ghk[23]=2;  Ghk[24]=0;  Ghk[25]=1;  Ghk[26]=2;
    Ghk[27]=0;  Ghk[28]=1;  Ghk[29]=2;  Ghk[30]=0;  Ghk[31]=1;  Ghk[32]=2;  Ghk[33]=0;  Ghk[34]=1;  Ghk[35]=2;  Ghk[36]=0;  Ghk[37]=1;  Ghk[38]=2;  Ghk[39]=0;  Ghk[40]=1;
    Ghk[41]=2;  Ghk[42]=0;  Ghk[43]=1;  Ghk[44]=2;  Ghk[45]=0;  Ghk[46]=1;  Ghk[47]=2;  Ghk[48]=0;  Ghk[49]=1;  Ghk[50]=2;  Ghk[51]=0;  Ghk[52]=1;  Ghk[53]=2;
    Ghk[54]=0;  Ghk[55]=1;  Ghk[56]=2;  Ghk[57]=0;  Ghk[58]=1;  Ghk[59]=2;  Ghk[60]=0;  Ghk[61]=1;  Ghk[62]=2;  Ghk[63]=0;  Ghk[64]=1;  Ghk[65]=2;  Ghk[66]=0;  Ghk[67]=1;
    Ghk[68]=2;  Ghk[69]=0;  Ghk[70]=1;  Ghk[71]=2;  Ghk[72]=0;  Ghk[73]=1;  Ghk[74]=2;  Ghk[75]=0;  Ghk[76]=1;  Ghk[77]=2;  Ghk[78]=0;  Ghk[79]=1;  Ghk[80]=2;
    return Ghk ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerVector Ghk(16) ;
    Ghk = IntegerVector::create (0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2) ;
    return Ghk ;
  }
  if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
    IntegerVector Ghk(36) ;
    Ghk[0]=0;   Ghk[1]=2;   Ghk[2]=0;   Ghk[3]=2;   Ghk[4]=0;   Ghk[5]=2;   Ghk[6]=0;   Ghk[7]=2;   Ghk[8]=0;   Ghk[9]=2;   Ghk[10]=0;  Ghk[11]=2;
    Ghk[12]=0;  Ghk[13]=2;  Ghk[14]=0;  Ghk[15]=2;  Ghk[16]=0;  Ghk[17]=2;  Ghk[18]=0;  Ghk[19]=2;  Ghk[20]=0;  Ghk[21]=2;  Ghk[22]=0;  Ghk[23]=2;
    Ghk[24]=0;  Ghk[25]=2;  Ghk[26]=0;  Ghk[27]=2;  Ghk[28]=0;  Ghk[29]=2;  Ghk[30]=0;  Ghk[31]=2;  Ghk[32]=0;  Ghk[33]=2;  Ghk[34]=0;  Ghk[35]=2;
    return Ghk ;
  }
  return 0;
}

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerMatrix trueGenotypes(int gender1, int gender2, int chrom){
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerMatrix emGeno(4,2) ;
    emGeno(0,0)=0; emGeno(0,1)=0;
    emGeno(1,0)=0; emGeno(1,1)=2;
    emGeno(2,0)=2; emGeno(2,1)=0;
    emGeno(3,0)=2; emGeno(3,1)=2;
    return emGeno ;
  }
  if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
    IntegerMatrix emGeno(6,2) ;
    emGeno(0,0)=0; emGeno(0,1)=0;
    emGeno(1,0)=0; emGeno(1,1)=2;
    emGeno(2,0)=1; emGeno(2,1)=0;
    emGeno(3,0)=1; emGeno(3,1)=2;
    emGeno(4,0)=2; emGeno(4,1)=0;
    emGeno(5,0)=2; emGeno(5,1)=2;
    return emGeno ;
  }
  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    IntegerMatrix emGeno(9,2) ;
    emGeno(0,0)=0; emGeno(0,1)=0;
    emGeno(1,0)=0; emGeno(1,1)=1;
    emGeno(2,0)=0; emGeno(2,1)=2;
    emGeno(3,0)=1; emGeno(3,1)=0;
    emGeno(4,0)=1; emGeno(4,1)=1;
    emGeno(5,0)=1; emGeno(5,1)=2;
    emGeno(6,0)=2; emGeno(6,1)=0;
    emGeno(7,0)=2; emGeno(7,1)=1;
    emGeno(8,0)=2; emGeno(8,1)=2;
    return emGeno ;
  }
  return 0;
}


//----- conditional emission probabilities -----

// Calculating the conditionaly emission probability for a single SNP
// @param ibd The IBD state of SNP j
// @param chrom The chromosome
// @param gender1 The gender of sample 1
// @param gender2 The gender of sample 2
// @param currentGeno The genotype of the current SNP
// @param conditionGeno The genotype of the condition SNP
// @param pAlleleFreq The population allele frequency
// @param pAlleleFreq The haplotype frequency
// @param error The genotype error rate
double conditionEmissionProb(const int ibd, const int chrom, const int gender1, const int gender2, IntegerVector currentGeno, IntegerVector conditionGeno, double pAlleleFreq, NumericVector hapltypeFreq, double error){
  NumericVector jointProb;
  IntegerVector Gij, Gik, Ghj, Ghk;
  IntegerMatrix trueGeno;
  double numerator = 0, denominator = 0;


  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    jointProb = jointProbDD(hapltypeFreq, ibd) ;
    Gij = GijC(gender1, gender2, chrom) ;
    Gik = GikC(gender1, gender2, chrom) ;
    Ghj = GhjC(gender1, gender2, chrom) ;
    Ghk = GhkC(gender1, gender2, chrom) ;
    for(int i = 0; i < 81; ++i){
      numerator += jointProb[i]*genotypeErrorD(Gij[i], conditionGeno[0] ,error)*genotypeErrorD(Gik[i],conditionGeno[1],error)*genotypeErrorD(Ghj[i],currentGeno[0],error)*genotypeErrorD(Ghk[i],currentGeno[1],error) ;
    }
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    jointProb = jointProbHH(hapltypeFreq, ibd) ;
    Gij = GijC(gender1, gender2, chrom) ;
    Gik = GikC(gender1, gender2, chrom) ;
    Ghj = GhjC(gender1, gender2, chrom) ;
    Ghk = GhkC(gender1, gender2, chrom) ;
    for(int i = 0; i < 16; ++i) numerator += jointProb[i]*genotypeErrorH(Gij[i], conditionGeno[0] ,error)*genotypeErrorH(Gik[i],conditionGeno[1],error)*genotypeErrorH(Ghj[i],currentGeno[0],error)*genotypeErrorH(Ghk[i],currentGeno[1],error) ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 2){
    jointProb = jointProbHD(hapltypeFreq, ibd) ;
    Gij = GijC(gender1, gender2, chrom) ;
    Gik = GikC(gender1, gender2, chrom) ;
    Ghj = GhjC(gender1, gender2, chrom) ;
    Ghk = GhkC(gender1, gender2, chrom) ;
    for(int i = 0; i < 36; ++i) numerator += jointProb[i]*genotypeErrorD(Gij[i], conditionGeno[1] ,error)*genotypeErrorH(Gik[i],conditionGeno[0],error)*genotypeErrorD(Ghj[i],currentGeno[1],error)*genotypeErrorH(Ghk[i],currentGeno[0],error) ;
  }
  if(chrom == 23 && gender1 == 2 && gender2 == 1){
    jointProb = jointProbHD(hapltypeFreq, ibd) ;
    Gij = GijC(gender1, gender2, chrom) ;
    Gik = GikC(gender1, gender2, chrom) ;
    Ghj = GhjC(gender1, gender2, chrom) ;
    Ghk = GhkC(gender1, gender2, chrom) ;
    for(int i = 0; i < 36; ++i) numerator += jointProb[i]*genotypeErrorD(Gij[i], conditionGeno[0] ,error)*genotypeErrorH(Gik[i],conditionGeno[1],error)*genotypeErrorD(Ghj[i],currentGeno[0],error)*genotypeErrorH(Ghk[i],currentGeno[1],error) ;
  }

  if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
    trueGeno = trueGenotypes(gender1, gender2, chrom) ;
    for(int g = 0; g < 9; ++g) denominator += emissionProbDD(pAlleleFreq, trueGeno(g,0), trueGeno(g,1), ibd)*genotypeErrorD(trueGeno(g,0),conditionGeno[0],error)*genotypeErrorD(trueGeno(g,1),conditionGeno[1],error) ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    trueGeno = trueGenotypes(gender1, gender2, chrom) ;
    for(int g = 0; g < 4; ++g) denominator += emissionProbHH(pAlleleFreq, trueGeno(g,0), trueGeno(g,1), ibd)*genotypeErrorH(trueGeno(g,0),conditionGeno[0],error)*genotypeErrorH(trueGeno(g,1),conditionGeno[1],error) ;
  }
  if(chrom == 23 && gender1 == 1 && gender2 == 2){
    trueGeno = trueGenotypes(gender1, gender2, chrom) ;
    for(int g = 0; g < 6; ++g) denominator += emissionProbHD(pAlleleFreq, trueGeno(g,0), trueGeno(g,1), ibd, 2, 1)*genotypeErrorD(trueGeno(g,0),conditionGeno[1],error)*genotypeErrorH(trueGeno(g,1),conditionGeno[0],error) ;
  }
  if(chrom == 23 && gender1 == 2 && gender2 == 1){
    trueGeno = trueGenotypes(gender1, gender2, chrom) ;
    for(int g = 0; g < 6; ++g) denominator += emissionProbHD(pAlleleFreq, trueGeno(g,0), trueGeno(g,1), ibd, 2, 1)*genotypeErrorD(trueGeno(g,0),conditionGeno[0],error)*genotypeErrorH(trueGeno(g,1),conditionGeno[1],error) ;
  }
  return numerator/denominator ;
}


//----- emission probability error summation -----

// Calculating the emission probability sumation when missing genotype calls present
// @param alleleFreq_t The population allele frequency of SNP t
// @param genotypes1_t The genotype of sample 1 at SNP t
// @param genotypes2_t The genotype of sample 2 at SNP t
// @param error The genotype error rate
// @param gender1 The gender of sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP t
// @param ibd_j The IBD state
double emissionProbMissingGeno(double alleleFreq_t, int genotypes1_t, int genotypes2_t, double error, int gender1, int gender2, int chrom, int ibd_j){
  IntegerMatrix trueGenotype ;
  trueGenotype = trueGenotypes(gender1,gender2,chrom) ; // All combinations of genotypes
  double emprobsError ;
  int noGenotypes;
  noGenotypes = trueGenotype.nrow() ; // total number of combinations of genotpes
  emprobsError = 0 ;

  for(int g = 0; g<noGenotypes; ++g){

    if(genotypes1_t == -1 && genotypes2_t != -1){ // first genotype missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbDD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),miss,error) * genotypeErrorD(trueGenotype(g,1),genotypes2_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),miss,error) * genotypeErrorD(trueGenotype(g,0),genotypes2_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss){
          emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),genotypes2_t,error) * genotypeErrorD(trueGenotype(g,0),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHH(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),miss,error) * genotypeErrorH(trueGenotype(g,1),genotypes2_t,error) ;
        }
      }
    }

    if(genotypes1_t != -1 && genotypes2_t == -1){ // second genotype missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbDD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorD(trueGenotype(g,1),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss = 0; miss < 3; miss = miss){
          emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),genotypes1_t,error) * genotypeErrorD(trueGenotype(g,0),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),miss,error) * genotypeErrorD(trueGenotype(g,0),genotypes1_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHH(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH(trueGenotype(g,1),miss,error) ;
        }
      }
    }

    if(genotypes1_t == -1 && genotypes2_t == -1){ // both genotypes missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emprobsError += emissionProbDD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),miss1,error) * genotypeErrorD(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),miss1,error) * genotypeErrorD(trueGenotype(g,0),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),miss1,error) * genotypeErrorH(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emprobsError += emissionProbHH(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),miss1,error) * genotypeErrorH(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
    }

    if(genotypes1_t != -1 && genotypes2_t != -1){ // no genotypes missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        emprobsError += emissionProbDD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorD(trueGenotype(g,1),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),genotypes1_t,error) * genotypeErrorD(trueGenotype(g,0),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        emprobsError += emissionProbHD(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH(trueGenotype(g,1),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        emprobsError += emissionProbHH(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH(trueGenotype(g,1),genotypes2_t,error) ;
      }
    }
  }
  return emprobsError ;
}


//----- conditional emission probability error summation -----

// Calculating the conditional emission probability sumation when missing genotype calls present
// @param alleleFreq_t The population allele frequency of SNP t
// @param current The genotypes of both samples at SNP t
// @param condition The genotypes of both samples at the condition SNP
// @param haplotypeFreq the haplotype frequency
// @param error The genotype error rate
// @param gender1 The gender of sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP t
// @param ibd_j The IBD state
double conditionEmissionProbMissingGeno(double alleleFreq_t, IntegerVector current, IntegerVector condition, NumericVector haplotypeFreq, double error, int gender1, int gender2, int chrom, int ibd_j){
  double emprobsCondition = 0.0;

  // no genotypes missing
  if(current[0] != -1 && current[1] != -1 && condition[0] != -1 && condition[1] != -1) emprobsCondition = conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;

  // single genotype missing
  // 1.
  if(current[0] == -1 && current[1] != -1 && condition[0] != -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2)){
      for(int miss = 0; miss < 3; ++miss){
        current[0] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
    if(chrom == 23 && gender1 == 1){
      for(int miss = 0; miss < 3; miss = miss + 2){
        current[0] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
  }
  // 2.
  if(current[0] != -1 && current[1] == -1 && condition[0] != -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender2 == 2)){
      for(int miss = 0; miss < 3; ++miss){
        current[1] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
    if(chrom == 23 && gender2 == 1){
      for(int miss = 0; miss < 3; miss = miss + 2){
        current[1] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
  }
  // 3.
  if(current[0] != -1 && current[1] != -1 && condition[0] == -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2)){
      for(int miss = 0; miss < 3; ++miss){
        condition[0] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
    if(chrom == 23 && gender1 == 1){
      for(int miss = 0; miss < 3; miss = miss + 2){
        condition[0] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
  }
  // 4.
  if(current[0] != -1 && current[1] != -1 && condition[0] != -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender2 == 2)){
      for(int miss = 0; miss < 3; ++miss){
        condition[1] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
    if(chrom == 23 && gender2 == 1){
      for(int miss = 0; miss < 3; miss = miss + 2){
        condition[1] = miss ;
        emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
      }
    }
  }
  // two genotypes missing
  // 1.
  if(current[0] == -1 && current[1] == -1 && condition[0] != -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[0] = miss1 ;
          current[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[0] = miss1 ;
          current[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[0] = miss1 ;
          current[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[0] = miss1 ;
          current[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // 2.
  if(current[0] == -1 && current[1] != -1 && condition[0] == -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[0] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[0] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // 3.
  if(current[0] == -1 && current[1] != -1 && condition[0] != -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // 4.
  if(current[0] != -1 && current[1] == -1 && condition[0] == -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[1] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[1] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[1] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[1] = miss1 ;
          condition[0] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // 5.
  if(current[0] != -1 && current[1] == -1 && condition[0] != -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          current[1] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          current[1] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // 6.
  if(current[0] != -1 && current[1] != -1 && condition[0] == -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          condition[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          condition[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          condition[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          condition[0] = miss1 ;
          condition[1] = miss2 ;
          emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
        }
      }
    }
  }
  // three genotypes missing
  // 1.
  if(current[0] == -1 && current[1] == -1 && condition[0] == -1 && condition[1] != -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[0] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[0] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[0] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[0] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
  }
  // 2.
  if(current[0] == -1 && current[1] == -1 && condition[0] != -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            current[1] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
  }
  // 3.
  if(current[0] != -1 && current[1] == -1 && condition[0] == -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[1] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[1] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[1] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[1] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
  }
  // 4.
  if(current[0] == -1 && current[1] != -1 && condition[0] == -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            current[0] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            current[0] = miss1 ;
            condition[0] = miss2 ;
            condition[1] = miss3 ;
            emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
          }
        }
      }
    }
  }
  // four genotypes missing
  // 1.
  if(current[0] == -1 && current[1] == -1 && condition[0] == -1 && condition[1] == -1){
    if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            for(int miss4 = 0; miss4 < 3; ++miss4){
              current[0] = miss1 ;
              current[1] = miss2 ;
              condition[0] = miss3 ;
              condition[1] = miss4 ;
              emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
            }
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 2){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; ++miss2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            for(int miss4 = 0; miss4 < 3; ++miss4){
              current[0] = miss1 ;
              current[1] = miss2 ;
              condition[0] = miss3 ;
              condition[1] = miss4 ;
              emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
            }
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 2 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; ++miss1){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; ++miss3){
            for(int miss4 = 0; miss4 < 3; miss4 = miss4 + 2){
              current[0] = miss1 ;
              current[1] = miss2 ;
              condition[0] = miss3 ;
              condition[1] = miss4 ;
              emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
            }
          }
        }
      }
    }
    if(chrom == 23 && gender1 == 1 && gender2 == 1){
      for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
        for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
          for(int miss3 = 0; miss3 < 3; miss3 = miss3 + 2){
            for(int miss4 = 0; miss4 < 3; miss4 = miss4 + 2){
              current[0] = miss1 ;
              current[1] = miss2 ;
              condition[0] = miss3 ;
              condition[1] = miss4 ;
              emprobsCondition += conditionEmissionProb(ibd_j, chrom, gender1, gender2, current, condition, alleleFreq_t, haplotypeFreq, error) ;
            }
          }
        }
      }
    }
  }
  return emprobsCondition ;
}


//----- alpha -----

// Calculate alpha
// @param noStates Integer. The number of IBD states in the model
// @param piProb A numeric vector containing the initial state probabilities
// @param meiosis Integer. The number of meiosis separating the two samples
// @param NoSNPs Integer. The number of SNPs
// @param genotypes A integer martix containing the genotype calls for a pair of samples
// @param alleleFreq A numeric vector of population allele frequencies
// @param positionM A numeric vector of SNP genetic map positions in M
// @param error Numeric. The genotype error rate
// @param gender_1 Integer. The gender of sample 1
// @param gender_2 Integer. The gender of sample 2
// @param conditionSNPs Integer vector. Numeric IDs of condition SNPs
// @param haplotypeFreq Numeric marix. haplotype frequencies
NumericMatrix calculate_alpha_m2(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom, IntegerVector conditionSNPs, NumericMatrix haplotypeFreq){
  // Initialising all parameters:
  NumericVector scale(NoSNPs);
  NumericVector alphaA(NoSNPs);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix alphaHat(NoSNPs, noStates);
  NumericVector  haplotypeVector(4);
  IntegerVector currentGeno(2), conditionGeno(2);
  double alphaASum, alphaSum, emprobsError, emprobsCondition;

  // Initialisation
  alphaSum = 0;
  for(int i = 0; i<noStates; ++i){
    emprobsError = emissionProbMissingGeno(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    alpha(0,i) = piProb[i] * emprobsError;
    alphaSum += alpha(0,i);
  }
  scale[0] = roundDecimal( 1/alphaSum, 4);
  for(int i = 0; i<noStates; ++i){
    alphaHat(0,i) = roundDecimal( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(NoSNPs-1); ++t){
    alphaSum = 0;
    for(int j = 0; j<noStates; ++j){
      alphaASum = 0;
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          alphaA[i] = alphaHat(t,i) * transitionProbDD(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1) ){
          alphaA[i] = alphaHat(t,i) * transitionProbHD(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          alphaA[i] = alphaHat(t,i) * transitionProbHH(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        alphaASum += alphaA[i];
      }

      currentGeno[0] = genotypes(t+1,0); currentGeno[1] = genotypes(t+1,1);
      conditionGeno[0] = genotypes(conditionSNPs[t+1],0); conditionGeno[1] = genotypes(conditionSNPs[t+1],1);
      haplotypeVector[0] = haplotypeFreq(t+1,0); haplotypeVector[1] = haplotypeFreq(t+1,1); haplotypeVector[2] = haplotypeFreq(t+1,2); haplotypeVector[3] = haplotypeFreq(t+1,3);
      emprobsCondition = conditionEmissionProbMissingGeno(alleleFreq[t+1], currentGeno, conditionGeno, haplotypeVector, error, gender1, gender2,chrom, j);

      alpha(t+1,j) = alphaASum * emprobsCondition;
      alphaSum += alpha(t+1,j);
    }
    scale[t+1] = 1/alphaSum;
    for(int i = 0; i<noStates; ++i){
      alphaHat(t+1,i) = roundDecimal( alpha(t+1,i) * scale[t+1], 4);
    }
  }
  return alphaHat  ;
}


//----- scale -----

// Calculate scale
// @param noStates Integer. The number of IBD states in the model
// @param piProb A numeric vector containing the initial state probabilities
// @param meiosis Integer. The number of meiosis separating the two samples
// @param NoSNPs Integer. The number of SNPs
// @param genotypes A integer martix containing the genotype calls for a pair of samples
// @param alleleFreq A numeric vector of population allele frequencies
// @param positionM A numeric vector of SNP genetic map positions in M
// @param error Numeric. The genotype error rate
// @param gender_1 Integer. The gender of sample 1
// @param gender_2 Integer. The gender of sample 2
// @param conditionSNPs Integer vector. Numeric IDs of condition SNPs
// @param haplotypeFreq Numeric marix. haplotype frequencies
NumericVector calculate_scale_m2(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom, IntegerVector conditionSNPs, NumericMatrix haplotypeFreq){
  // Initialising all parameters:
  NumericVector scale(NoSNPs);
  NumericVector alphaA(NoSNPs);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix alphaHat(NoSNPs, noStates);
  NumericVector  haplotypeVector(4);
  IntegerVector currentGeno(2), conditionGeno(2);
  double alphaASum, alphaSum, emprobsError, emprobsCondition;
  emprobsCondition = 0.0;

  // Initialisation
  alphaSum = 0;
  for(int i = 0; i<noStates; ++i){
    emprobsError = emissionProbMissingGeno(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    alpha(0,i) = piProb[i] * emprobsError;
    alphaSum += alpha(0,i);
  }
  scale[0] = roundDecimal( 1/alphaSum, 4);
  for(int i = 0; i<noStates; ++i){
    alphaHat(0,i) = roundDecimal( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(NoSNPs-1); ++t){
    alphaSum = 0;
    for(int j = 0; j<noStates; ++j){
      alphaASum = 0;
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          alphaA[i] = alphaHat(t,i) * transitionProbDD(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1) ){
          alphaA[i] = alphaHat(t,i) * transitionProbHD(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          alphaA[i] = alphaHat(t,i) * transitionProbHH(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        alphaASum += alphaA[i];
      }

      currentGeno[0] = genotypes(t+1,0); currentGeno[1] = genotypes(t+1,1);
      conditionGeno[0] = genotypes(conditionSNPs[t+1],0); conditionGeno[1] = genotypes(conditionSNPs[t+1],1);
      haplotypeVector[0] = haplotypeFreq(t+1,0); haplotypeVector[1] = haplotypeFreq(t+1,1); haplotypeVector[2] = haplotypeFreq(t+1,2); haplotypeVector[3] = haplotypeFreq(t+1,3);
      emprobsCondition = conditionEmissionProbMissingGeno(alleleFreq[t+1], currentGeno, conditionGeno, haplotypeVector, error, gender1, gender2,chrom, j);
      alpha(t+1,j) = alphaASum * emprobsCondition;
      alphaSum += alpha(t+1,j);
    }
    scale[t+1] = 1/alphaSum;
    for(int i = 0; i<noStates; ++i){
      alphaHat(t+1,i) = roundDecimal( alpha(t+1,i) * scale[t+1], 4);
    }
  }
  return scale  ;
}


//----- beta -----

// Calculate beta
// @param noStates Integer. The number of IBD states in the model
// @param piProb A numeric vector containing the initial state probabilities
// @param meiosis Integer. The number of meiosis separating the two samples
// @param NoSNPs Integer. The number of SNPs
// @param genotypes A integer martix containing the genotype calls for a pair of samples
// @param alleleFreq A numeric vector of population allele frequencies
// @param positionM A numeric vector of SNP genetic map positions in M
// @param scale Numeric vector of scales used in alpha
// @param error Numeric. The genotype error rate
// @param gender_1 Integer. The gender of sample 1
// @param gender_2 Integer. The gender of sample 2
// @param conditionSNPs Integer vector. Numeric IDs of condition SNPs
// @param haplotypeFreq Numeric marix. haplotype frequencies
NumericMatrix calculate_beta_m2(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, NumericVector scale, double error, int gender1, int gender2, int chrom, IntegerVector conditionSNPs, NumericMatrix haplotypeFreq){
  // Initialising all parameters:
  NumericMatrix beta(NoSNPs, noStates);
  NumericMatrix betaHat(NoSNPs, noStates);
  IntegerVector currentGeno(2), conditionGeno(2);
  NumericVector  haplotypeVector(4);
  double betaSum, emprobsCondition;
  int conditionSNP;
  const int T = NoSNPs-1;

  // Initialisation
  for(int i = 0; i<noStates; ++i){
    beta(T,i) = 1;
    betaHat(T,i) = roundDecimal( beta(T,i) * scale[T], 4);
  }

  // Induction
  for(int t = (T-1); t>=0; t--){
    for(int i = 0; i<noStates; ++i){
      betaSum = 0;
      for(int j = 0; j<noStates; ++j){
        if(t == 0){
          emprobsCondition = emissionProbMissingGeno(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, j);
        }
        if(t != 0){
          currentGeno[0] = genotypes(t+1,0); currentGeno[1] = genotypes(t+1,1);
          conditionGeno[0] = genotypes(conditionSNPs[t+1],0); conditionGeno[1] = genotypes(conditionSNPs[t+1],1);
          haplotypeVector[0] = haplotypeFreq(t+1,0); haplotypeVector[1] = haplotypeFreq(t+1,1); haplotypeVector[2] = haplotypeFreq(t+1,2); haplotypeVector[3] = haplotypeFreq(t+1,3);
          emprobsCondition = conditionEmissionProbMissingGeno(alleleFreq[t+1], currentGeno, conditionGeno, haplotypeVector, error, gender1, gender2, chrom, j);
        }
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          beta(t,j) = transitionProbDD(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsCondition * betaHat(t+1,j) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          beta(t,j) = transitionProbHH(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsCondition * betaHat(t+1,j) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
          beta(t,j) = transitionProbHD(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsCondition * betaHat(t+1,j) ;
        }
        betaSum += beta(t,j);
      }
      betaHat(t,i) = betaSum * scale[t];
    }
    double betaHatSum ;
    if(noStates == 2){
      betaHatSum = betaHat(t,0) + betaHat(t,1) ;
      betaHat(t,0) = betaHat(t,0)/betaHatSum ;
      betaHat(t,1) = betaHat(t,1)/betaHatSum ;
    }
    if(noStates == 3){
      betaHatSum = betaHat(t,0) + betaHat(t,1) + betaHat(t,2) ;
      betaHat(t,0) = betaHat(t,0)/betaHatSum ;
      betaHat(t,1) = betaHat(t,1)/betaHatSum ;
      betaHat(t,2) = betaHat(t,2)/betaHatSum ;
    }

  }
  return betaHat;
}


//----- Viterbi -----

// Calculate viterbi
// @param noStates Integer. The number of IBD states in the model
// @param piProb A numeric vector containing the initial state probabilities
// @param meiosis Integer. The number of meiosis separating the two samples
// @param NoSNPs Integer. The number of SNPs
// @param genotypes A integer martix containing the genotype calls for a pair of samples
// @param alleleFreq A numeric vector of population allele frequencies
// @param positionM A numeric vector of SNP genetic map positions in M
// @param error Numeric. The genotype error rate
// @param gender_1 Integer. The gender of sample 1
// @param gender_2 Integer. The gender of sample 2
// @param conditionSNPs Integer vector. Numeric IDs of condition SNPs
// @param haplotypeFreq Numeric marix. haplotype frequencies
// [[Rcpp::export]]
IntegerVector calculate_viterbi_m2(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom, IntegerVector conditionSNPs, NumericMatrix haplotypeFreq){
  // Initialising all parameters:
  IntegerVector qStar(NoSNPs);
  NumericVector deltaA(NoSNPs);
  NumericMatrix delta(NoSNPs, noStates);
  NumericMatrix psi(NoSNPs, noStates);
  IntegerVector currentGeno(2), conditionGeno(2);
  NumericVector haplotypeVector(4);
  double emprobsError, emprobsCondition;
  double tProb;
  const int T = NoSNPs-1;

  // Initialisation
  for(int i = 0; i<noStates; ++i){
    psi(0,i) = 0;
    emprobsError = emissionProbMissingGeno(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    tProb = log(piProb[i]);
    if(tProb < -100000) tProb = log(0.00001); // stop log of -infinity
    delta(0,i) = tProb + log(emprobsError) ;
  }

  // Recursion
  for(int t = 1; t<NoSNPs; ++t){
    for(int j = 0; j<noStates; ++j){
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          tProb = log(transitionProbDD(piProb[0],piProb[1],piProb[2],meiosis,positionM[t]-positionM[t-1],j,i)) ;
          if(tProb < -100000) tProb = log(0.00001);
          deltaA[i] = delta(t-1,i) + tProb ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          tProb = log(transitionProbHH(piProb[0],meiosis,positionM[t]-positionM[t-1],j,i)) ;
          if(tProb < -100000) tProb = log(0.00001);
          deltaA[i] = delta(t-1,i) + tProb ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
          tProb = log(transitionProbHD(piProb[0],meiosis,positionM[t]-positionM[t-1],j,i)) ;
          if(tProb < -100000) tProb = log(0.00001);
          deltaA[i] = delta(t-1,i) + tProb ;
        }
      }
      currentGeno[0] = genotypes(t,0); currentGeno[1] = genotypes(t,1);
      conditionGeno[0] = genotypes(conditionSNPs[t],0); conditionGeno[1] = genotypes(conditionSNPs[t],1);
      haplotypeVector[0] = haplotypeFreq(t,0); haplotypeVector[1] = haplotypeFreq(t,1); haplotypeVector[2] = haplotypeFreq(t,2); haplotypeVector[3] = haplotypeFreq(t,3);
      emprobsCondition = conditionEmissionProbMissingGeno(alleleFreq[t], currentGeno, conditionGeno, haplotypeVector, error, gender1, gender2,chrom, j);

      double maxNum = deltaA[0];
      int maxArg = 0;
      for(int k = 0; k<noStates; ++k){
        if(deltaA[k]>maxNum){
          maxNum = deltaA[k];
          maxArg = k; // could be problematic if there are more than 1 maximum
        }
      }
      delta(t,j) = maxNum + log(emprobsCondition) ;
      psi(t,j) = maxArg ;
    }
  }

  // Termination
  qStar[T] = 0;
  double probStar = delta(T,0);
  for(int k = 0; k<noStates; ++k){
    if(delta(T,k)>probStar){
      qStar[T] = k; // could be problematic if there are more than 1 maximum
      probStar = delta(T,k);
    }
  }

  // Path backtracking
  for(int t = (T-1); t>=0; t--){
    qStar[t] = psi(t+1,qStar[t+1]); // Not sure if this will work without qstar being a constant
  }
  return qStar;
}


//----- gamma -----

// Calculate gamma
// @param noStates Integer. The number of IBD states in the model
// @param piProb A numeric vector containing the initial state probabilities
// @param meiosis Integer. The number of meiosis separating the two samples
// @param NoSNPs Integer. The number of SNPs
// @param genotypes A integer martix containing the genotype calls for a pair of samples
// @param alleleFreq A numeric vector of population allele frequencies
// @param positionM A numeric vector of SNP genetic map positions in M
// @param error Numeric. The genotype error rate
// @param gender_1 Integer. The gender of sample 1
// @param gender_2 Integer. The gender of sample 2
// @param conditionSNPs Integer vector. Numeric IDs of condition SNPs
// @param haplotypeFreq Numeric marix. haplotype frequencies
// [[Rcpp::export]]
NumericMatrix calculate_gamma_m2(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom, IntegerVector conditionSNPs, NumericMatrix haplotypeFreq){
  NumericMatrix gamma(NoSNPs, noStates);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix beta(NoSNPs, noStates);
  NumericVector scale(NoSNPs);
  NumericVector alphaBeta(NoSNPs);
  double alphaBetaSum;

  alpha = calculate_alpha_m2(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, error, gender1, gender2, chrom, conditionSNPs, haplotypeFreq);
  scale = calculate_scale_m2(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, error, gender1, gender2, chrom, conditionSNPs, haplotypeFreq);
  beta  = calculate_beta_m2(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, scale, error, gender1, gender2, chrom, conditionSNPs, haplotypeFreq);

  for(int t = 0; t<NoSNPs; ++t){
    alphaBetaSum = 0;
    for(int i = 0; i<noStates; ++i){
      alphaBeta[i] = alpha(t,i) * beta(t,i);
      alphaBetaSum += alphaBeta[i];
    }
    for(int i = 0; i<noStates; ++i){
      gamma(t,i) = roundDecimal( (alpha(t,i) * beta(t,i))/alphaBetaSum, 3) ;
    }
  }
  return gamma;
}


