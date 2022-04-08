#' Internal Function
#'
#' IBD Summary Table
#'
#' Extracts the SNP with the largest LD with the current SNP
#' and returns this as the condition SNP
#' @param genotypes A data frame containing genotype information for a set of SNPs.
#' See \code{\link{getGenotypes}} for more details.
#' @param snp.ld A data frame generated from PLINK containing information on LD between pairs
#' of SNPs (A and B). This data frame contains exactly 7 columns of informations:
#' \enumerate{
#' \item Chromosome of SNP A (type \code{"numeric"} or \code{"integer"})
#' \item Base-pair position of SNP A (type \code{"numeric"} or \code{"integer"})
#' \item SNP A identifier (type \code{"character"})
#' \item Chromosome of SNP B (type \code{"numeric"} or \code{"integer"})
#' \item Base-pair position of SNP B (type \code{"numeric"} or \code{"integer"})
#' \item SNP B identifier (type \code{"character"})
#' \item R-sequared LD statistic between SNP A and SNP B (type \code{"numeric"})
#' }
#' where each row contains the LD information for a single pair of SNPs. The data frame should contain the header
#' \code{colnames(snp.ld)=c(CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2)}
#' @param maximum.ld.r2 A numeric value denoting the maximum linkage disequilibrium R2 value allowed
#' between pairs of SNPs.
#' @return A data frame containing columns
#' \enumerate{
#' \item Chromosome (type \code{numeric} or \code{integer})
#' \item SNP identifier (type \code{character})
#' \item Genetic map position (type \code{numeric} or \code{integer})
#' \item Base-pair position (type \code{numeric} or \code{integer})
#' \item Numeric SNP identifier (type \code{numeric} or \code{integer})
#' \item Numeric condition SNP identifier (type \code{numeric} or \code{integer})
#' }
#' with columns headed \code{chr, snp_id, pos_M, pos_bp, marker.id} and \code{cond.snp.id}.
getConditionSNP <- function(genotypes, snp.ld, maximum.ld.r2) {
  # create initial set of condition SNPs to be updated
  condition.snps <- data.frame(genotypes[,c("chr","snp_id","pos_M","pos_bp")],
                               marker.id = rep(0,nrow(genotypes)), cond.snp.id = rep(0,nrow(genotypes)))
  for (i in unique(condition.snps[,"chr"])) {
    condition.snps.chr <- condition.snps[condition.snps[,"chr"] == i,]
    condition.snps[condition.snps[,"chr"] == i,"marker.id"] <- 1:nrow(condition.snps.chr)
    condition.snps[condition.snps[,"chr"] == i,"cond.snp.id"][2:nrow(condition.snps.chr)] <- 1:(nrow(condition.snps.chr)-1)
  }

  # subset LD by SNPs in genotypes and R2 less than threshold
  snp.ld.1 <- snp.ld[snp.ld[,"SNP_A"] %in% genotypes[,"snp_id"] & snp.ld[,"SNP_B"] %in% genotypes[,"snp_id"] & snp.ld[,"R2"] <= maximum.ld.r2,]

  # for each unique SNP, find the SNP with the largest
  if (nrow(snp.ld.1) > 0) {
    condition.snps.0 <- NULL
    for (i in unique(snp.ld.1[,"SNP_B"])) {
      snp.ld.i <- snp.ld.1[snp.ld.1[,"SNP_B"] == i,]
      cond.snp <- snp.ld.i[snp.ld.i[,"R2"] == max(snp.ld.i[,"R2"]),"SNP_A"][1]
      cond.snp.id <- condition.snps[condition.snps[,"snp_id"] == cond.snp, "marker.id"]
      condition.snps.0 <- rbind(condition.snps.0,cbind(i, cond.snp.id))
    }
    colnames(condition.snps.0) <- c("snp_id","cond.snp.id")

    # merge with initial dataset
    condition.snps.1 <- merge(condition.snps,condition.snps.0,by="snp_id",all.x=TRUE)
    condition.snps.1[!is.na(condition.snps.1[,"cond.snp.id.y"]),"cond.snp.id.x"] <- as.character(condition.snps.1[!is.na(condition.snps.1[,"cond.snp.id.y"]),"cond.snp.id.y"])
    condition.snps.1 <- condition.snps.1[order(condition.snps.1[,"chr"],condition.snps.1[,"pos_bp"]),]

    # select columns of interest
    condition.snps.2 <- condition.snps.1[,c("chr","snp_id","pos_M","pos_bp","marker.id","cond.snp.id.x")]
    colnames(condition.snps.2) <- c("chr","snp_id","pos_M","pos_bp","marker.id","cond.snp.id")

  } else {
    condition.snps.2 <- condition.snps
  }
  condition.snps.2[,"marker.id"] <- as.numeric(condition.snps.2[,"marker.id"])
  condition.snps.2[,"cond.snp.id"] <- as.numeric(condition.snps.2[,"cond.snp.id"])

  return(condition.snps.2)
}


#' Calculates haplotype frequencies between two SNPs
#' using cpp scripts.
#' @param genotypes A data frame containing genotype information for a set of SNPs.
#' See \code{\link{getGenotypes}} for more details.
#' @param condition.snps A dataframe containging information on SNPs to condition on.
#' See \code{\link{getConditionSNP}} for more details.
#' @param genders A numeric vector of the genders of the individuals in the genotypes
#' data frame. 1 = male, 2 = female.
#' @return see return value from \code{\link{getGenotypes}} for more details.
getHaplotypeFreq <- function(genotypes, condition.snps, genders) {
  # create matrix for results
  condition.matrix <- matrix(nrow=nrow(genotypes),ncol=6)

  # for each chromosome, get the haplotype frequencies
  start.row <- 1
  for (i in unique(condition.snps[,"chr"])) {
    genotypes.chr <- as.matrix(genotypes[genotypes[,"chr"] == i,6:ncol(genotypes)])
    condition.snps.chr <- as.matrix(condition.snps[condition.snps[,"chr"] == i,c("marker.id","cond.snp.id")])
    condition.matrix[start.row:(start.row + nrow(genotypes.chr) - 1),] <- calculateHaplotypeFreq(genotypes.chr, condition.snps.chr, genders, i, polyroot)
    start.row <- start.row + nrow(genotypes.chr)
  }
  # add column names and check for NAs
  condition.matrix.0 <- cbind(genotypes[,c("chr","snp_id","pos_M","pos_bp","freq")], condition.matrix)
  colnames(condition.matrix.0) <- c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")
  if (any(is.na(condition.matrix.0))) stop("An error has occurred when determining haplotype frequencies. Use 'model=1'")
  return(condition.matrix.0)
}




