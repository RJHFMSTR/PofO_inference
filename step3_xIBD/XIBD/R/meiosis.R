# Internal Function
#
# Estimation of Meiosis
#
# calculateMeiosis() estimates the number of meiosis separating a pair of isolates given the global IBD pop_allele_freqs estimates.
# This method is described in Purcell et al (2007).
# @param omega.0 A numeric value between 0 and 1 representing the pop_allele_freqs of sharing 0 alleles IBD. The sum of omega.0, omega.1 and omega.2 should equal 1.
# @param omega.1 A numeric value between 0 and 1 representing the pop_allele_freqs of sharing 1 alleles IBD.
# @param omega.2 A numeric value between 0 and 1 representing the pop_allele_freqs of sharing 2 alleles IBD.
# @return The number of meiosis separating the pair of isoaltes.
calculateMeiosis <- function(omega.0,omega.1,omega.2){
  if(omega.1==0&omega.2==0){
    m <- 20
  }
  if(omega.1!=0&omega.2==0){
    z1 <- omega.1
    z2 <- 0
    m <- 1 - (log(z1)/log(2))
  }
  if(omega.2!=0){
    coef1 <- omega.2
    coef2 <- -(omega.1 + 2*omega.2)
    coef3 <- 1
    roots <- Re(polyroot(c(coef1,coef2,coef3)))

    z1 <- roots[1]
    z2 <- omega.2/z1
    #z2 <- omega.1/z1
    m1 <- 1 - (log(z1)/log(2))
    m2 <- 1 - (log(z2)/log(2))
    #m <- min(m1,m2)[1]
    m <- m1 + m2
  }
  return(m)
}


# IBDparameters() calculates IBD probabilities then estimates the number of meiosis for an isolate pair.
# @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
# @param pop_allele_freqs A numeric vector of population allele frequencies for each SNP.
# @param gender_1 An integer denoting the MOI value of isoalte 1.
# @param gender_2 An integer denoting the MOI value of isoalte 2.
# @return A vector of 4 values representing the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD respectively.
IBDparameters <- function(genotypes, pop_allele_freqs, gender_1, gender_2, chromosome){
  if (chromosome == 23 & gender_1 == 1 & gender_2 == 1) {
    b <- bVectorHH(genotypes)
    A <- AmatrixHH(pop_allele_freqs,genotypes)
    x <- solve(t(A), b)
    x <- round(x,4)
    x[3] <- 0
    if (x[1] > 1) {
      x[1] = 1
      x[2] = 0
    }
  }
  if (chromosome == 23 & ((gender_1 == 1 & gender_2 == 2) | (gender_1 == 2 & gender_2 == 1))) {
    b <- bVectorHD(genotypes)
    A <- AmatrixHD(pop_allele_freqs,genotypes)
    x <- solve(t(A), b)
    x <- round(x,4)
    x[3] <- 0
    if (x[1] > 1) {
      x[1] = 1
      x[2] = 0
    }
  }
  if (chromosome != 23 | (gender_1 == 2 & gender_2 == 2)) {
    b <- bVectorDD(genotypes)
    A <- AmatrixDD(pop_allele_freqs,genotypes)
    x <- solve(t(A), b)
    if (x[1] > 1) {
      x[1] = 1
      x[2] = 0
      x[3] = 0
    }
    if (x[3] < 0) {
      x[3] = 0
      s = x[1] + x[2]
      x[1] = x[1]/s
      x[2] = x[2]/s
    }
    pi. <- x[2]/2 + x[3]
    if (pi.^2 <= x[3]) {
      x[1] = (1 - pi.)^2
      x[2] = 2*pi.*(1-pi.)
      x[3] = pi.^2
    }
    x <- round(x,4)
  }
  meiosis <- calculateMeiosis(x[1], x[2], x[3])
  ibd_estimates <- cbind(meiosis, x[1], x[2], x[3])
  colnames(ibd_estimates) <- c("m","z0","z1","z2")
  return(ibd_estimates)
}
