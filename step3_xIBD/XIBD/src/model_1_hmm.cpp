#include <Rcpp.h>
using namespace Rcpp;

// Internal Function
//
// HMM Model 1
//
// All the functions needed from the HMM, including emission probabilities
// transition probabilities, alpha, beta, gamma and viterbi.
// Model is based on Purcell et al. 2007


//----- round decimal -----

// Round digits to specified decimal places
// @param number A number to round
// @param number The number of digits to round to
double roundDecimal1(double number, int digits){
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
double emissionProbHH1(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  double alt_allele_freq = 1 - pop_allele_freq;
  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal1( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal1( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal1( pow(alt_allele_freq,2), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal1( pop_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal1( alt_allele_freq, 6);

  return 0;
}


// The emission probabilities for 1 haploid chromosome and 1 diploid chromosome
// @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
// @param genotype_1 The genotype for isolate 1 from the pair for SNP i
// @param genotype_2 The genotype for isolate 2 from the pair for SNP i
// @param ibd The IBD state
// @param male_column The haploid isolate from the pair. Either 1 or 2
// @param female_column The diploid isolate from the pair. Either 1 or 2
double emissionProbHD1(double pop_allele_freq, int genotype_1, int genotype_2, int ibd, int male_column, int female_column) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  int geno_male, geno_female;

  if(male_column == 1) { geno_male = genotype_1; geno_female = genotype_2; }
  if(male_column == 2) { geno_male = genotype_2; geno_female = genotype_1; }

  if(geno_female == 0 && geno_male == 0 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,3), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 0 && ibd == 0) return roundDecimal1( 2*pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 0) return roundDecimal1( 2*pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 0) return roundDecimal1( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 2 && ibd == 0) return roundDecimal1( pow(alt_allele_freq,3), 6);

  if(geno_female == 0 && geno_male == 0 && ibd == 1) return roundDecimal1( pow(pop_allele_freq,2), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 1) return 0;
  if(geno_female == 1 && geno_male == 0 && ibd == 1) return roundDecimal1( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 1) return roundDecimal1( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 1) return 0;
  if(geno_female == 2 && geno_male == 2 && ibd == 1) return roundDecimal1( pow(alt_allele_freq,2), 6);

  return 0;
}


// The emission probabilities for 2 diploid chromosomes
// @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
// @param genotype_1 The genotype for isolate 1 from the pair for SNP i
// @param genotype_2 The genotype for isolate 2 from the pair for SNP i
// @param ibd The IBD state
double emissionProbDD1(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,4), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 0) return roundDecimal1( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 0) return roundDecimal1( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 0) return roundDecimal1( 4*pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 0) return roundDecimal1( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal1( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 0) return roundDecimal1( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal1( pow(alt_allele_freq,4), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal1( pow(pop_allele_freq,3), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 1) return roundDecimal1( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 1) return roundDecimal1( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 1) return roundDecimal1( (pow(pop_allele_freq,2)*alt_allele_freq + pop_allele_freq*pow(alt_allele_freq,2)), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 1) return roundDecimal1( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 1) return roundDecimal1( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal1( pow(alt_allele_freq,3), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 2) return roundDecimal1( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 2) return roundDecimal1( 2*pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 2) return roundDecimal1( pow(alt_allele_freq,2), 6);

  return 0;
}


//----- transition probabilities -----

// The transition probabilities for 2 haploid chromosomes
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbHH1(double omega_0, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal1( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal1( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal1( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal1( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}


// The transition probabilities for 1 haploid and 1 diploid chromosome
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbHD1(double omega_0, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal1( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal1( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal1( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal1( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}


// The transition probabilities for 2 diploid chromosomes
// @param omega_0 The probability of sharing 0 alleles IBD
// @param meiosis The number of meiosis separating the two isoaltes
// @param dist_cM The genetic map distance (cM) between SNP i and SNP j
// @param ibd_current The IBD state of SNP j
// @param ibd_previous The IBD state of SNP i
double transitionProbDD1(double omega_0, double omega_1, double omega_2, int meiosis, double dist_M, int ibd_current, int ibd_previous) {
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  double T02 = (exp(-alpha * omega_1 * dist_M) * omega_2)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_0 * omega_1)/(omega_1 - 1) + omega_2;
  double T20 = (exp(-alpha * omega_1 * dist_M) * omega_0)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_2 * omega_1)/(omega_1 - 1) + omega_0;

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal1( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T02), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal1( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 0 && ibd_current == 2) return roundDecimal1( T02, 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal1( ((1 - exp(-alpha * dist_M)) * omega_0), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal1( ((1 - exp(-alpha * dist_M)) * omega_1 + exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 2) return roundDecimal1( ((1 - exp(-alpha * dist_M)) * omega_2), 6);
  if(ibd_previous == 2 && ibd_current == 0) return roundDecimal1( T20, 6);
  if(ibd_previous == 2 && ibd_current == 1) return roundDecimal1( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 2 && ibd_current == 2) return roundDecimal1( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T20), 6);

  return 0;
}


//----- genotype error probabilities -----

// The genotyping error probability for 1 haploid chromosome
// @param truth The true genotype
// @param observed The observed genotype
// @param error The genotype error rate
double genotypeErrorH1(int truth, int observed, double error){
  if(truth == observed) return 1 - error;
  else return error;
}


// The genotyping error probability for 1 diploid chromosome
// @param truth The true genotype
// @param observed The observed genotype
// @param error The genotype error rate
double genotypeErrorD1(int truth, int observed, double error){
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


//----- true genotypes -----

// Matrices of all possible genotype combinations between pairs, given genders
// @param gender1 The gender if sample 1
// @param gender2 The gender of sample 2
// @param chrom The chromosome of SNP i
IntegerMatrix trueGenotypes1(int gender1, int gender2, int chrom){
  // 2 male X chromosomes
  if(chrom == 23 && gender1 == 1 && gender2 == 1){
    IntegerMatrix emGeno(4,2) ;
    emGeno(0,0)=0; emGeno(0,1)=0;
    emGeno(1,0)=0; emGeno(1,1)=2;
    emGeno(2,0)=2; emGeno(2,1)=0;
    emGeno(3,0)=2; emGeno(3,1)=2;
    return emGeno ;
  }
  // 1 male and 1 female X chromosome
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
  // 2 female X chromosomes or 2 autosomes
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
double emissionProbMissingGeno1(double alleleFreq_t, int genotypes1_t, int genotypes2_t, double error, int gender1, int gender2, int chrom, int ibd_j){
  IntegerMatrix trueGenotype ;
  trueGenotype = trueGenotypes1(gender1,gender2,chrom) ; // All combinations of genotypes
  double emprobsError ;
  int noGenotypes, miss, miss1, miss2 ;
  noGenotypes = trueGenotype.nrow() ; // total number of combinations of genotpes
  emprobsError = 0 ;

  for(int g = 0; g<noGenotypes; ++g){

    if(genotypes1_t == -1 && genotypes2_t != -1){ // first genotype missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbDD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD1(trueGenotype(g,0),miss,error) * genotypeErrorD1(trueGenotype(g,1),genotypes2_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH1(trueGenotype(g,1),miss,error) * genotypeErrorD1(trueGenotype(g,0),genotypes2_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD1(trueGenotype(g,0),miss,error) * genotypeErrorH1(trueGenotype(g,1),genotypes2_t,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHH1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH1(trueGenotype(g,0),miss,error) * genotypeErrorH1(trueGenotype(g,1),genotypes2_t,error) ;
        }
      }
    }

    if(genotypes1_t != -1 && genotypes2_t == -1){ // second genotype missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbDD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorD1(trueGenotype(g,1),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH1(trueGenotype(g,1),genotypes1_t,error) * genotypeErrorD1(trueGenotype(g,0),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss = 0; miss < 3; ++miss){
          emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH1(trueGenotype(g,1),miss,error) ;
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emprobsError += emissionProbHH1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH1(trueGenotype(g,1),miss,error) ;
        }
      }
    }

    if(genotypes1_t == -1 && genotypes2_t == -1){ // both genotypes missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emprobsError += emissionProbDD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD1(trueGenotype(g,0),miss1,error) * genotypeErrorD1(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH1(trueGenotype(g,1),miss1,error) * genotypeErrorD1(trueGenotype(g,0),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD1(trueGenotype(g,0),miss1,error) * genotypeErrorH1(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emprobsError += emissionProbHH1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH1(trueGenotype(g,0),miss1,error) * genotypeErrorH1(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
    }

    if(genotypes1_t != -1 && genotypes2_t != -1){ // no genotypes missing
      if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
        emprobsError += emissionProbDD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorD1(trueGenotype(g,1),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 2){
        emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH1(trueGenotype(g,1),genotypes1_t,error) * genotypeErrorD1(trueGenotype(g,0),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 2 && gender2 == 1){
        emprobsError += emissionProbHD1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH1(trueGenotype(g,1),genotypes2_t,error) ;
      }
      if(chrom == 23 && gender1 == 1 && gender2 == 1){
        emprobsError += emissionProbHH1(alleleFreq_t,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH1(trueGenotype(g,0),genotypes1_t,error) * genotypeErrorH1(trueGenotype(g,1),genotypes2_t,error) ;
      }
    }
  }
  return emprobsError ;
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
// @param chrom The chromosome
NumericMatrix calculate_alpha_m1(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom){
  // Initialising all parameters:
  NumericVector scale(NoSNPs);
  NumericVector alphaA(NoSNPs);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix alphaHat(NoSNPs, noStates);
  double alphaASum, alphaSum, emprobsError;

  // Initialisation
  alphaSum = 0;
  for(int i = 0; i<noStates; ++i){
    emprobsError = emissionProbMissingGeno1(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    alpha(0,i) = piProb[i] * emprobsError;
    alphaSum += alpha(0,i);
  }
  scale[0] = roundDecimal1( 1/alphaSum, 4);
  for(int i = 0; i<noStates; ++i){
    alphaHat(0,i) = roundDecimal1( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(NoSNPs-1); ++t){
    alphaSum = 0;
    for(int j = 0; j<noStates; ++j){
      alphaASum = 0;
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          alphaA[i] = alphaHat(t,i) * transitionProbDD1(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1) ){
          alphaA[i] = alphaHat(t,i) * transitionProbHD1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          alphaA[i] = alphaHat(t,i) * transitionProbHH1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        alphaASum += alphaA[i];
      }
      emprobsError = emissionProbMissingGeno1(alleleFreq[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender1, gender2, chrom, j) ;
      alpha(t+1,j) = alphaASum * emprobsError;
      alphaSum += alpha(t+1,j);
    }
    scale[t+1] = roundDecimal1( 1/alphaSum, 4);
    for(int i = 0; i<noStates; ++i){
      alphaHat(t+1,i) = roundDecimal1( alpha(t+1,i) * scale[t+1], 4);
    }
  }
  return alphaHat;
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
// @param chrom The chromosome
NumericVector calculate_scale_m1(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom){
  // Initialising all parameters:
  NumericVector scale(NoSNPs);
  NumericVector alphaA(NoSNPs);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix alphaHat(NoSNPs, noStates);
  double alphaASum, alphaSum, emprobsError;

  // Initialisation
  alphaSum = 0;
  for(int i = 0; i<noStates; ++i){
    emprobsError = emissionProbMissingGeno1(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    alpha(0,i) = piProb[i] * emprobsError;
    alphaSum += alpha(0,i);
  }
  scale[0] = roundDecimal1( 1/alphaSum, 4);
  for(int i = 0; i<noStates; ++i){
    alphaHat(0,i) = roundDecimal1( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(NoSNPs-1); ++t){
    alphaSum = 0;
    for(int j = 0; j<noStates; ++j){
      alphaASum = 0;
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          alphaA[i] = alphaHat(t,i) * transitionProbDD1(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1) ){
          alphaA[i] = alphaHat(t,i) * transitionProbHD1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          alphaA[i] = alphaHat(t,i) * transitionProbHH1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) ;
        }
        alphaASum += alphaA[i];
      }
      emprobsError = emissionProbMissingGeno1(alleleFreq[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender1, gender2, chrom, j) ;
      alpha(t+1,j) = alphaASum * emprobsError;
      alphaSum += alpha(t+1,j);
    }
    scale[t+1] = roundDecimal1( 1/alphaSum, 4);
    for(int i = 0; i<noStates; ++i){
      alphaHat(t+1,i) = roundDecimal1( alpha(t+1,i) * scale[t+1], 4);
    }
  }
  return scale;
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
// @param chrom The chromosome
NumericMatrix calculate_beta_m1(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, NumericVector scale, double error, int gender1, int gender2, int chrom){
  // Initialising all parameters:
  NumericMatrix beta(NoSNPs, noStates);
  NumericMatrix betaHat(NoSNPs, noStates);
  IntegerMatrix trueGenotype;
  int noGenotypes;
  double betaSum, emprobsError;
  const int T = NoSNPs-1;

  trueGenotype = trueGenotypes1(gender1,gender2,chrom) ; // All combinations of genotypes
  noGenotypes  = trueGenotype.nrow() ; // total number of combinations of genotpes

  // Initialisation
  for(int i = 0; i<noStates; ++i){
    beta(T,i)    = 1;
    betaHat(T,i) = roundDecimal1( beta(T,i) * scale[T], 4);
  }

  // Induction
  for(int t = (T-1); t>=0; t--){
    for(int i = 0; i<noStates; ++i){
      betaSum = 0;
      for(int j = 0; j<noStates; ++j){
        emprobsError = emissionProbMissingGeno1(alleleFreq[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender1, gender2, chrom, j) ;
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          beta(t,j) = transitionProbDD1(piProb[0],piProb[1],piProb[2],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsError * betaHat(t+1,j) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          beta(t,j) = transitionProbHH1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsError * betaHat(t+1,j) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
          beta(t,j) = transitionProbHD1(piProb[0],meiosis,positionM[t+1]-positionM[t],j,i) * emprobsError * betaHat(t+1,j) ;
        }
        betaSum += beta(t,j);
      }
      betaHat(t,i) = roundDecimal1( betaSum * scale[t], 4);
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
// @param chrom The chromosome
// [[Rcpp::export]]
IntegerVector calculate_viterbi_m1(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom){
  // Initialising all parameters:
    IntegerVector qStar(NoSNPs);
  NumericVector deltaA(NoSNPs);
  NumericMatrix delta(NoSNPs, noStates);
  NumericMatrix psi(NoSNPs, noStates);
  IntegerMatrix trueGenotype;
  int noGenotypes;
  double emprobsError;
  const int T = NoSNPs-1;

  // Initialisation
  for(int i = 0; i<noStates; ++i){
    psi(0,i) = 0;
    emprobsError = emissionProbMissingGeno1(alleleFreq[0], genotypes(0,0), genotypes(0,1), error, gender1, gender2, chrom, i) ;
    delta(0,i) = log(piProb[i]) + log(emprobsError) ;
  }

  // Recursion
  for(int t = 1; t<NoSNPs; ++t){
    for(int j = 0; j<noStates; ++j){
      for(int i = 0; i<noStates; ++i){
        if(chrom != 23 || (chrom == 23 && gender1 == 2 && gender2 == 2)){
          deltaA[i] = delta(t-1,i) + log(transitionProbDD1(piProb[0],piProb[1],piProb[2],meiosis,positionM[t]-positionM[t-1],j,i)) ;
        }
        if(chrom == 23 && gender1 == 1 && gender2 == 1){
          deltaA[i] = delta(t-1,i) + log(transitionProbHH1(piProb[0],meiosis,positionM[t]-positionM[t-1],j,i)) ;
        }
        if((chrom == 23 && gender1 == 1 && gender2 == 2) || (chrom == 23 && gender1 == 2 && gender2 == 1)){
          deltaA[i] = delta(t-1,i) + log(transitionProbHD1(piProb[0],meiosis,positionM[t]-positionM[t-1],j,i)) ;
        }
      }
      emprobsError = emissionProbMissingGeno1(alleleFreq[t], genotypes(t,0), genotypes(t,1), error, gender1, gender2, chrom, j) ;
      double maxNum = deltaA[0];
      int maxArg = 0;
      for(int k = 0; k<noStates; ++k){
        if(deltaA[k]>maxNum){
          maxNum = deltaA[k];
          maxArg = k; // could be problematic if there are more than 1 maximum
        }
      }
      delta(t,j) = maxNum + log(emprobsError) ;
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
// @param chrom The chromosome
// [[Rcpp::export]]
NumericMatrix calculate_gamma_m1(const int noStates, NumericVector piProb, int meiosis, const int NoSNPs, IntegerMatrix genotypes, NumericVector alleleFreq, NumericVector positionM, double error, int gender1, int gender2, int chrom){
  NumericMatrix gamma(NoSNPs, noStates);
  NumericMatrix alpha(NoSNPs, noStates);
  NumericMatrix beta(NoSNPs, noStates);
  NumericVector scale(NoSNPs);
  NumericVector alphaBeta(NoSNPs);
  double alphaBetaSum;

  alpha = calculate_alpha_m1(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, error, gender1, gender2, chrom);
  scale = calculate_scale_m1(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, error, gender1, gender2, chrom);
  beta  = calculate_beta_m1(noStates, piProb, meiosis, NoSNPs, genotypes, alleleFreq, positionM, scale, error, gender1, gender2, chrom);

  for(int t = 0; t<NoSNPs; ++t){
    alphaBetaSum = 0;
    for(int i = 0; i<noStates; ++i){
      alphaBeta[i] = alpha(t,i) * beta(t,i);
      alphaBetaSum += alphaBeta[i];
    }
    for(int i = 0; i<noStates; ++i){
      gamma(t,i) = roundDecimal1( (alpha(t,i) * beta(t,i))/alphaBetaSum, 3) ;
    }
  }
  return gamma;
}



