#' to prevent notes
globalVariables("pair.i")
#' XIBD Parameter Estimation
#'
#' Estimate the number of meioses and the probabilities of sharing 0, 1 and 2 alleles IBD between pairs.
#'
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param number.cores the number of cores used for parallel execution.
#' @return A named list containing \code{autosome_parameters} and \code{X_chromosome_parameters}.
#' Both \code{autosome_parameters} and \code{X_chromosome_parameters} are data frames with columns:
#' \enumerate{
#' \item Family 1 ID (type \code{"character"})
#' \item Individual 1 ID (type \code{"character"})
#' \item Family 2 ID (type \code{"character"})
#' \item Individual 2 ID (type \code{"character"})
#' \item The number of meiosis (type \code{"numeric"} or \code{"integer"})
#' \item Probability of sharing 0 alleles IBD (type \code{"numeric"})
#' \item Probability of sharing 1 allele IBD (type \code{"numeric"})
#' \item Probability of sharing 2 alleles IBD (type \code{"numeric"})
#' }
#' headed \code{fid1, iid1, fid2, iid2, m, ibd0, ibd1} and \code{ibd2}, respectively.
#' Each row describes parameters for a unique pair of samples.
#' \code{autosome_parameters} will be \code{NULL} if chromosomes 1-22 are all excluded
#' from \code{ped.genotypes} and \code{X_chromosome_parameters} will be \code{NULL} if
#' chromosome 23 is excluded from \code{ped.genotypes}.
#' @importFrom foreach "%dopar%"
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile
#' @export
#' @examples
#' # estimate parameters
#' my_parameters <- getIBDparameters(ped.genotypes = example_genotypes,
#'                                   number.cores = 1)
#'
#' str(my_parameters)
getIBDparameters <- function(ped.genotypes, number.cores = 1){

  # check input ped and genotypes
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  stopifnot(c("pedigree","genotypes","model") %in% names(ped.genotypes))
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  model     <- ped.genotypes[["model"]]

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) < 8 & nrow(genotypes) <= 1)
    stop ("ped.genotypes has incorrect format")
  colnames(genotypes)[1:5] <- c("chr", "snp_id", "pos_M","pos_bp", "freq")

  # check numeric input parameters
  stopifnot(is.numeric(number.cores))

  # check haploid status
  # stopifnot(is.logical(haploid))

  # create new isolate IDs from PED FIDs and IIDs
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  number.pairs  <- 1:nrow(isolate.pairs)

  # calculate quantiles from the number of pairs, for progress bar only
  pair.quantiles   <- unique(round(quantile(number.pairs,probs=seq(0,0.9999,0.01))))
  number.quantiles <- length(pair.quantiles)

  # create progress bar
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)

  # define number of cores
  doParallel::registerDoParallel(cores=number.cores)

  # for each subgroup of pairs belonging to the quantiles, get parameters (for progress bar)
  start <- 1
  ibd.estimates <- list()
  for (quantile.group in 1:number.quantiles) {

    # assign pairs to a subgroup based on which quantile they're in (for progress bar)
    if (number.quantiles == nrow(isolate.pairs))
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      pair.group <- start:(start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles))
      pair.group <- start:nrow(isolate.pairs)

    # get IBD parameters for subgroups of pairs
    ibd.estimates.0 <- foreach::foreach(pair.i=pair.group, .combine='mergeLists1') %dopar% {
      fid.1    <- as.character(isolate.pairs[pair.i,1])
      iid.1    <- as.character(isolate.pairs[pair.i,2])
      fid.2    <- as.character(isolate.pairs[pair.i,3])
      iid.2    <- as.character(isolate.pairs[pair.i,4])
      gender.1 <- pedigree[pedigree[,"fid"] == fid.1 & pedigree[,"iid"] == iid.1,"sex"]
      gender.2 <- pedigree[pedigree[,"fid"] == fid.2 & pedigree[,"iid"] == iid.2,"sex"]
      pair.genotypes   <- cbind(genotypes[,paste(fid.1,iid.1,sep="/")], genotypes[,paste(fid.2,iid.2,sep="/")])
      pop.allele.freqs <- genotypes[,"freq"]

      #if (!haploid) {
      # autosome parameters
      ibd.estimates.1 <- NULL
      if (any(genotypes[,"chr"] != 23)) {
        pair.genotypes.1 <- pair.genotypes[genotypes[,"chr"] != 23,]
        pop.allele.freqs.1 <- pop.allele.freqs[genotypes[,"chr"] != 23]
        ibd.estimates.1 <- cbind(fid.1, iid.1, fid.2, iid.2, IBDparameters(pair.genotypes.1, pop.allele.freqs.1, gender.1, gender.2, 1))
        colnames(ibd.estimates.1) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
      }
      # X chromosome parameters
      ibd.estimates.23 <- NULL
      if (any(genotypes[,"chr"] == 23)) {
        pair.genotypes.23 <- pair.genotypes[genotypes[,"chr"] == 23,]
        pop.allele.freqs.23 <- pop.allele.freqs[genotypes[,"chr"] == 23]
        ibd.estimates.23 <- cbind(fid.1, iid.1, fid.2, iid.2, IBDparameters(pair.genotypes.23, pop.allele.freqs.23, gender.1, gender.2, 23))
        colnames(ibd.estimates.23) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
      }
      # create list of parameters
      ibd.estimates <- list(ibd.estimates.1, ibd.estimates.23)

      #} #else {
      # haploid parameters
      #ibd.estimates <- cbind(fid.1, iid.1, fid.2, iid.2, IBDparameters(pair.genotypes, pop.allele.freqs, 1, 1, 23))
      #colnames(ibd.estimates) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
      #ibd.estimates <- list(ibd.estimates)
      #}
      ibd.estimates
    }

    # merge parameter estimates for this subgroup with the global set
    ibd.estimates <- mergeLists1(ibd.estimates, ibd.estimates.0)

    # assign new start pair for next subgroup
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      start <- start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group]

    # update progress bar
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)

  # format parameter estimates
  if (any(genotypes[,"chr"] != 23) & !any(genotypes[,"chr"] == 23)) {
    if (length(ibd.estimates) == 1)
      ibd.estimates <- list(ibd.estimates[[1]], NULL)
    ibd.estimates[[1]] <- data.frame(ibd.estimates[[1]])
    for (i in 1:4)
      ibd.estimates[[1]][,i] <- as.character(ibd.estimates[[1]][,i])
    for (i in 5:8)
      ibd.estimates[[1]][,i] <- as.numeric(as.character(ibd.estimates[[1]][,i]))
  }
  if (!any(genotypes[,"chr"] != 23) & any(genotypes[,"chr"] == 23)) {
    if (length(ibd.estimates) == 1)
      ibd.estimates <- list(NULL, ibd.estimates[[1]])
    ibd.estimates[[2]] <- data.frame(ibd.estimates[[2]])
    for (i in 1:4)
      ibd.estimates[[2]][,i] <- as.character(ibd.estimates[[2]][,i])
    for (i in 5:8)
      ibd.estimates[[2]][,i] <- as.numeric(as.character(ibd.estimates[[2]][,i]))
  }
  if (any(genotypes[,"chr"] != 23) & any(genotypes[,"chr"] == 23)){
    for (j in 1:length(ibd.estimates)) {
      ibd.estimates[[j]] <- data.frame(ibd.estimates[[j]])
      for (i in 1:4)
        ibd.estimates[[j]][,i] <- as.character(ibd.estimates[[j]][,i])
      for (i in 5:8)
        ibd.estimates[[j]][,i] <- as.numeric(as.character(ibd.estimates[[j]][,i]))
    }
  }
  names(ibd.estimates) <- c("autosome_parameters","X_chromosome_parameters")

  return(ibd.estimates)
}
