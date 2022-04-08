#' to prevent notes
globalVariables("pair.i")
#' XIBDs IBD Segment Detection
#'
#' Detects genomic regions shared IBD between pairs.
#'
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' Note the family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param parameters a data frame containing meioses and IBD probability estimates for all pairwise combinations of samples.
#' See \code{Value} description in \code{\link{getIBDparameters}} for more details.
#' @param model an integer of either 1 or 2 denoting which of the two models should be run.
#' \enumerate{
#' \item \code{model=1} is based on the HMM implemented in PLINK (Purcell et al., 2007) which assumes the SNPs
#' are in linkage equilibrium (LE). This often requires thinning of markers prior to use.
#' \item \code{model=2} is based on the HMM implemented in RELATE (Albrechtsen et al., 2009) which allows SNPs
#' to be in LD and implicitly accounts for the LD through conditional emission probabilities where
#' the current genotype probability is conditioned on the genotype of a single previous SNP.
#' }
#' The default is \code{model=NULL} and the model listed in \code{ped.genotypes} is used. NOTE: \code{model=2} can
#' only be used if \code{model=2} is in \code{ped.genotypes}.
#' @param chromosomes a numeric vector containing a subset of chromosomes to perform IBD analysis on. The
#' default is \code{chromosomes=NULL} and IBD analysis will be performed on all chromosomes in \code{ped.genotypes}.
#' @param number.cores the number of cores used for parallel execution.
#' @param minimum.snps the minimum number of SNPs in an IBD segment for it to be reported. The default value is 20 SNPs.
#' @param minimum.length.bp the minimum length of a reported IBD segment. The default value is 50,000 bp.
#' @param error the genotyping error rate. The default value is 0.001.
#' @param posterior a logical value indicating whether posterior probabilities for each pairwise analysis should be returned.
#' The posterior probability is calculated for each SNP as posteriorPr(IBD=1)/2 + posteriorPr(IBD=2) using the
#' forward and backward variables (Rabiner, 1989).
#' A data frame containing probabilities for each SNP and each pairwise analysis is returned.
#' This data frame can be very large when there are many SNPs and many pairwise analyses, and the run-time of
#' \code{getIBDsegments()} will increase. The default is \code{posterior=FALSE}; \code{posterior=TRUE} is not recommended for large datasets.
#' @return A named list of 1 object when \code{posterior=FALSE} and 2 objects when \code{posterior=TRUE}.
#' The first object in the list, \code{ibd_segments}, is a data frame with information:
#' \enumerate{
#' \item Family 1 ID (type \code{"character"})
#' \item Individual 1 ID (type \code{"character"})
#' \item Family 2 ID (type \code{"character"})
#' \item Individual 2 ID (type \code{"character"})
#' \item Chromosome (type \code{"numeric"} or \code{integer})
#' \item SNP identifier (type \code{"character"})
#' \item Start SNP (type \code{"character"})
#' \item End SNP (type \code{"character"})
#' \item Start position bp (type \code{"numeric"} or \code{integer})
#' \item End position bp (type \code{"numeric"} or \code{integer})
#' \item Start position M (type \code{"numeric"})
#' \item End position M (type \code{"numeric"})
#' \item Number of SNPs (type \code{"numeric"} or \code{integer})
#' \item Length bp (type \code{"numeric"} or \code{integer})
#' \item Length M (type \code{"numeric"})
#' \item IBD status (1 = one allele shared IBD, 2 = two alleles shared IBD) (type \code{"numeric"} or \code{integer})
#' }
#' where each row is a unique IBD segment for a pair of individuals. The data frame is headed
#' \code{fid1, iid1, fid2, iid2, chr, start.snp, end.snp, start.position.bp, end.position.bp, start.position.M, end.position.M,
#' number.snps, length.bp, length.M, ibd.status}. The second object (returned when \code{posterior=TRUE}), \code{posterior_probabilities}, is a
#' data frame with the first four columns
#' \enumerate{
#' \item Chromosome (type \code{"numeric"} or \code{integer})
#' \item SNP identifier (type \code{"character"})
#' \item Genetic map distance (type \code{"numeric"})
#' \item Base-pair positions (type \code{"numeric"} or \code{integer})
#' }
#' and columns 5 onwards are the posterior probabilities for each pair with pair identifier headers.
#' Rows correspond to SNPs.
#' @importFrom foreach "%dopar%"
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile
#' @export
#' @examples
#' \dontrun{
#' # infer IBD
#' my_ibd <- getIBDsegments(ped.genotypes = example_genotypes,
#'                          parameters = example_parameters,
#'                          model = NULL,
#'                          chromosomes = NULL,
#'                          number.cores = 1,
#'                          minimum.snps = 20,
#'                          minimum.length.bp = 50000,
#'                          error = 0.001,
#'                          posterior = FALSE)
#'
#' str(my_ibd)
#' }
getIBDsegments <- function(ped.genotypes, parameters, model = NULL, chromosomes = NULL, number.cores = 1,
                           minimum.snps = 20, minimum.length.bp = 50000, error = 0.001, posterior = FALSE){

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  stopifnot(c("pedigree","genotypes","model") %in% names(ped.genotypes))
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  model_0   <- ped.genotypes[["model"]]
  if (is.null(model))
    model <- model_0
  stopifnot(is.numeric(model))
  stopifnot(model == 1 | model == 2)

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (model == 1){
    if(ncol(genotypes) < 8 & nrow(genotypes) <= 1)
      stop ("ped.genotypes has incorrect format")
    if(!all(colnames(genotypes)[1:5] %in% c("chr", "snp_id", "pos_M","pos_bp", "freq")))
      stop ("ped.genotypes has incorrect format")
  }
  if (model == 2) {
    if(ncol(genotypes) < 14 & nrow(genotypes) <= 1)
      stop ("ped.genotypes has incorrect format")
    if(!all(colnames(genotypes)[1:11] %in% c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")))
      stop ("ped.genotypes has incorrect format")
  }

  # check paramters file is a dataframe with correct fields
  stopifnot(is.list(parameters) & length(parameters) == 2)
  for (i in 1:length(parameters)) {
    if(!is.null(parameters[[i]])) {
      if (ncol(parameters[[i]]) != 8)
        stop ("parameters has incorrect format")
      colnames(parameters[[i]]) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
    }
    if(!all(c("autosome_parameters","X_chromosome_parameters") %in% names(parameters)))
      stop("'autosome_parameters' and 'X_chromosome_parameters' are missing from 'parameters'")
  }
  parameters.a <- parameters[["autosome_parameters"]]
  parameters.x <- parameters[["X_chromosome_parameters"]]

  # check chromosomes
  if (!is.null(chromosomes)) {
    stopifnot(is.numeric(chromosomes))
    if(!all(chromosomes %in% genotypes[,"chr"]))
      stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% genotypes[,"chr"])]," not in 'ped.genotypes'\n")))
  } else
    chromosomes <- as.numeric(unique(as.character(genotypes[,"chr"])))
  chromosomes <- chromosomes[order(chromosomes)]
  if (any(chromosomes < 23) & is.null(parameters.a))
    stop("no parameters estimated for autosomes")
  if (any(chromosomes == 23) & is.null(parameters.x))
    stop("no parameters estimated for the X chromosome")

  # check numeric input parameters
  stopifnot(is.numeric(number.cores))
  stopifnot(is.numeric(minimum.snps))
  stopifnot(is.numeric(minimum.length.bp))
  stopifnot(is.numeric(error))
  stopifnot(is.logical(posterior))

  # create new isolate IDs from PED FIDs and IIDs
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  number.pairs  <- 1:nrow(isolate.pairs)

  # calculate quantiles from the number of pairs, for progress bar only
  pair.quantiles <- unique(round(quantile(number.pairs,probs=seq(0,0.9999,0.01))))
  number.quantiles   <- length(pair.quantiles);

  # create progress bar
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)

  # define number of cores
  doParallel::registerDoParallel(cores=number.cores)

  # for each subgroup of pairs belonging to the quantiles, get IBD segments (for progress bar)
  start <- 1
  ibd.segments <- list()
  for (quantile.group in 1:number.quantiles) {

    # assign pairs to a subgroup based on which quantile they're in (for progress bar)
    if (number.quantiles == nrow(isolate.pairs))
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      pair.group <- start:(start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles))
      pair.group <- start:nrow(isolate.pairs)

    # get IBD segments for subgroups of pairs
    ibd.segments.0 <- foreach::foreach(pair.i=pair.group, .combine='mergeLists2') %dopar% {

      # select pair
      fid.1    <- as.character(isolate.pairs[pair.i,1])
      iid.1    <- as.character(isolate.pairs[pair.i,2])
      fid.2    <- as.character(isolate.pairs[pair.i,3])
      iid.2    <- as.character(isolate.pairs[pair.i,4])
      gender.1 <- pedigree[pedigree[,"fid"] == fid.1 & pedigree[,"iid"] == iid.1,"sex"]
      gender.2 <- pedigree[pedigree[,"fid"] == fid.2 & pedigree[,"iid"] == iid.2,"sex"]

      # for each chromosome, perform IBD analysis
      ibd.table.2 <- NULL
      posterior.prob <- NULL
      for (chrom in chromosomes) {
        # define number of states in model based on genders
        if (chrom != 23 | (gender.1 == 2 & gender.2 == 2)) {
          number.states = 3
        } else
          number.states = 2

        # get model paramters
        if (chrom != 23 & !is.null(parameters.a)) {
          meiosis      <- as.numeric(parameters.a[parameters.a[,"fid1"] == fid.1 & parameters.a[,"fid2"] == fid.2 & parameters.a[,"iid1"] == iid.1 & parameters.a[,"iid2"] == iid.2,"m"])
          initial.prob <- as.numeric(parameters.a[parameters.a[,"fid1"] == fid.1 & parameters.a[,"fid2"] == fid.2 & parameters.a[,"iid1"] == iid.1 & parameters.a[,"iid2"] == iid.2,c("ibd0","ibd1","ibd2")])
        }
        if(chrom == 23 & !is.null(parameters.x)) {
          meiosis      <- as.numeric(parameters.x[parameters.x[,"fid1"] == fid.1 & parameters.x[,"fid2"] == fid.2 & parameters.x[,"iid1"] == iid.1 & parameters.x[,"iid2"] == iid.2,"m"])
          initial.prob <- as.numeric(parameters.x[parameters.x[,"fid1"] == fid.1 & parameters.x[,"fid2"] == fid.2 & parameters.x[,"iid1"] == iid.1 & parameters.x[,"iid2"] == iid.2,c("ibd0","ibd1","ibd2")])
        }

        # change initial probabilities to allow switiching between states.. re-think?? FIXME!!
        if (initial.prob[1] == 1) {
          initial.prob[1] <- 0.999
          initial.prob[2] <- 0.001
        }
        if (initial.prob[2] == 1) {
          initial.prob[1] <- 0.001
          initial.prob[2] <- 0.999
        }
        if (sum(initial.prob) != 1) {
          initial.prob <- c(initial.prob[1]/sum(initial.prob), initial.prob[2]/sum(initial.prob), initial.prob[3]/sum(initial.prob))
        }

        pair.genotypes <- cbind(genotypes[genotypes[,"chr"] == chrom,paste(fid.1,iid.1,sep="/")], genotypes[genotypes[,"chr"] == chrom,paste(fid.2,iid.2,sep="/")])
        positions.m  <- genotypes[genotypes[,"chr"] == chrom,"pos_M"]
        positions.bp <- genotypes[genotypes[,"chr"] == chrom,"pos_bp"]
        chromosome   <- as.character(genotypes[genotypes[,"chr"] == chrom,"chr"])
        markers      <- as.character(genotypes[genotypes[,"chr"] == chrom,"snp_id"])
        number.snps  <- length(positions.m)
        gamma        <- NULL

        # simple HMM - model 1
        if(model == 1){
          pop.allele.freqs <- genotypes[genotypes[,"chr"] == chrom,"freq"]
          viterbi <- calculate_viterbi_m1(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2, chrom)
          if (posterior)
            gamma <- calculate_gamma_m1(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2, chrom)
        }
        # conditional HMM - model 2
        if(model == 2){
          pop.allele.freqs <- genotypes[genotypes[,"chr"] == chrom,"freq_condition_snp"]
          condition.snps   <- genotypes[genotypes[,"chr"] == chrom,"condition_snp"] - 1
          haplotype.freqs  <- as.matrix(genotypes[genotypes[,"chr"] == chrom,c("pba","pbA","pBa","pBA")],ncol=4,nrow=number.snps)
          viterbi <- calculate_viterbi_m2(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2, chrom, condition.snps, haplotype.freqs)
          if (posterior)
            gamma <- calculate_gamma_m2(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2, chrom, condition.snps, haplotype.freqs)
        }

        # format posterior probabilities
        if (posterior) {
          if(number.states == 2) gamma <- cbind(gamma, rep(0,dim(gamma)[1]))
          posterior.prob <- c(posterior.prob, (gamma[,2]/2 + gamma[,3]))
        }

        # get IBD summary table
        ibd.results <- cbind(fid.1, iid.1, fid.2, iid.2, 1:number.snps, chromosome, markers, positions.m, positions.bp, viterbi)
        colnames(ibd.results) <- c("fid1","iid1","fid2","iid2","markerNo","chr","marker","pos.m","pos.bp","viterbi")
        ibd.table.1 <- IBDTable(ibd.results)

        # remove IBD segments with less than minimum.snps and less than minimum.length.bp
        # FIXME! can take a while to rbind
        if(length(ibd.table.1) != 0){
          ibd.table.2 <- rbind(ibd.table.2, ibd.table.1[as.numeric(ibd.table.1[,"number.snps"]) >= minimum.snps & as.numeric(ibd.table.1[,"length.bp"]) >= minimum.length.bp,])
        }
      }

      list(ibd.table.2, posterior.prob)
    }
    ibd.segments <- mergeLists2(ibd.segments, ibd.segments.0)

    # assign new start pair for next subgroup
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      start <- start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group]

    # update progress bar
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)

  # format IBD segments
  if(length(ibd.segments) > 0) {
    if (length(ibd.segments) == 2) {
      if (!is.null(ibd.segments[[1]])) {
        rownames(ibd.segments[[1]]) <- NULL
        ibd.segments[[1]] <- data.frame(ibd.segments[[1]])
        for(i in 1:7)
          ibd.segments[[1]][,i] <- as.character(ibd.segments[[1]][,i])
        for(i in 8:15)
          ibd.segments[[1]][,i] <- as.numeric(as.character(ibd.segments[[1]][,i]))

        number.pairs.ibd <- length(unique(paste(ibd.segments[[1]][,1],ibd.segments[[1]][,2],ibd.segments[[1]][,3],ibd.segments[[1]][,4])))
        cat(paste(number.pairs.ibd,"pairs inferred IBD\n"))
        cat(paste(nrow(ibd.segments[[1]]),"IBD segments detected\n"))
      }
      if (!is.null(ibd.segments[[2]])) {
        genotypes.chr <- NULL
        for (chrom in chromosomes) {
          genotypes.chr <- rbind(genotypes.chr, genotypes[genotypes[,"chr"] == chrom,c("chr","snp_id","pos_M","pos_bp")])
        }
        colnames(ibd.segments[[2]]) <- paste0(isolate.pairs[,1],".",isolate.pairs[,2],"/",isolate.pairs[,3],".",isolate.pairs[,4])
        ibd.segments[[2]] <- cbind(genotypes.chr, ibd.segments[[2]])
      }
      names(ibd.segments) <- c("ibd_segments","posterior_probabilities")
    }

    if (length(ibd.segments) == 1 & !posterior) {
      rownames(ibd.segments[[1]]) <- NULL
      ibd.segments[[1]] <- data.frame(ibd.segments[[1]])
      for(i in 1:7)
        ibd.segments[[1]][,i] <- as.character(ibd.segments[[1]][,i])
      for(i in 8:15)
        ibd.segments[[1]][,i] <- as.numeric(as.character(ibd.segments[[1]][,i]))
      names(ibd.segments)[1] <- "ibd_segments"

      number.pairs.ibd <- length(unique(paste(ibd.segments[[1]][,1],ibd.segments[[1]][,2],ibd.segments[[1]][,3],ibd.segments[[1]][,4])))
      cat(paste(number.pairs.ibd,"pairs inferred IBD\n"))
      cat(paste(nrow(ibd.segments[[1]]),"IBD segments detected\n"))
    }

    if (length(ibd.segments) == 1 & posterior) {
      genotypes.chr <- NULL
      for (chrom in chromosomes) {
        genotypes.chr <- rbind(genotypes.chr, genotypes[genotypes[,"chr"] == chrom,c("chr","snp_id","pos_M","pos_bp")])
      }
      colnames(ibd.segments[[1]]) <- paste0(isolate.pairs[,1],".",isolate.pairs[,2],"/",isolate.pairs[,3],".",isolate.pairs[,4])
      ibd.segments[[1]] <- cbind(genotypes.chr, ibd.segments[[1]])
      ibd.segments <- list(NULL, ibd.segments[[1]])
      names(ibd.segments) <- c("ibd_segments","posterior_probabilities")
    }
  } else
    cat("No IBD segments detected")

  return(ibd.segments)
}


