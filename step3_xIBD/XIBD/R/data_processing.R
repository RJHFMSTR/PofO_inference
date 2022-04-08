#' XIBD Pre-Analysis Data Processing
#'
#' Perform pre-analysis data processing of PLINK formatted unphased haplotype data,
#' including removal of SNPs and samples with high proportions of missing data, SNPs with low minor
#' allele frequencies and SNPs in high linkage disequilibrium (LD, based on R^2 if \code{model=2}).
#' Also, calculate population allele frequencies as well as haplotype frequencies between pairs of
#' SNPs (based on R^2 if \code{model=2}).
#'
#' @param ped.map a list with 2 objects:
#' \enumerate{
#' \item a data frame which contains the PLINK PED information.
#' The first six columns of this data frame are:
#' \enumerate{
#' \item Family ID (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item Individual ID (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item Paternal ID (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item Maternal ID (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item Gender (1 = male, 2 = female)
#' \item Phenotype (1 = unaffected, 2 = affected, 0 = unknown)
#' }
#' where each row describes a single sample The IDs are alphanumeric: the combination of family and individual ID should uniquely identify a sample. The phenotype column is not used
#' in XIBD analyses however it is required for completeness of a standard pedigree.
#' Columns 7 onwards are the sample haplotypes where the A and B alleles are coded as 1 and 2 respectively and missing data is coded as 0.
#' All SNPs (whether haploid or not) must have two alleles specified and each allele should be in a separate column. For example,
#' the alleles in columns 7 and 8 correspond to the unphased haplotypes of SNP 1 in the map file.
#' For haploid chromosomes, haplotypes should be specified as homozygous. Either both alleles should be missing (i.e. 0)
#' or neither. No header row should be given.
#' \item a data frame which contains the PLINK MAP information. This data frame contains exactly four columns of information:
#' \enumerate{
#' \item Chromosome (\code{"numeric"} or \code{"integer"})
#' \item SNP identifier (type \code{"character"})
#' \item Genetic map distance (centi morgans cM, or morgans M - default) (type \code{"numeric"})
#' \item Base-pair position (type \code{"numeric"} or \code{"integer"})
#' }
#' where each row describes a single marker. Genetic map distance and base-pair positions are expected to be positive values. The MAP file must
#' be ordered by increasing chromosomes and positions. SNP identifiers can contain any characters expect spaces or tabs; also you should avoid
#' * symbols in the names. The MAP file must contain as many markers as are in the PED file. No header row should be given.
#' }
#' @param reference.ped.map a list containing reference data used to calculate population allele
#' frequencies and haplotype frequencies, in the same format as \code{ped.map}.
#' The default value is \code{reference.ped.map=NULL} and XIBD will calculate the population allele
#' frequencies and haplotype frequencies from \code{ped.map}. This is not recommended for small datasets or
#' datasets of mixed populations. Genetic map positions and base-pair positions in this dataset are used
#' in the analysis. HapMap phase 2 and 3 PED and MAP data (hg19/build 37) for the 11 HapMap populations
#' can be downloaded from \url{http://bioinf.wehi.edu.au/software/XIBD}.
#' @param snp.ld optional for \code{model=1}; compulsory for \code{model=2}.
#' A data frame generated from PLINK containing information on LD between pairs of SNPs (A and B). This
#' data frame contains exactly 7 columns of information:
#' \enumerate{
#' \item Chromosome of SNP A (type \code{"numeric"} or \code{"integer"})
#' \item Base-pair position of SNP A (type \code{"numeric"} or \code{"integer"})
#' \item SNP A identifier (type \code{"character"})
#' \item Chromosome of SNP B (type \code{"numeric"} or \code{"integer"})
#' \item Base-pair position of SNP B (type \code{"numeric"} or \code{"integer"})
#' \item SNP B identifier (type \code{"character"})
#' \item R-squared LD statistic between SNP A and SNP B (type \code{"numeric"})
#' }
#' where each row contains the LD information for a single pair of SNPs. The data frame should contain the header
#' \code{CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B} and \code{R2}.
#' HapMap phase 2 and 3 PED and MAP data (hg19/build 37) for the 11 HapMap populations
#' can be downloaded from \url{http://bioinf.wehi.edu.au/software/XIBD}.
#' Alternatively, an LD file can be created using using PLINK (\url{http://www.cog-genomics.org/plink2}).
#' @param model an integer of either 1 or 2 denoting which of the two models should be run.
#' \enumerate{
#' \item \code{model=1} is based on the HMM implemented in PLINK (Purcell et al., 2007) which assumes the SNPs
#' are in linkage equilibrium (LE). This often requires thinning of markers prior to use.
#' \item \code{model=2} is based on the HMM implemented in RELATE (Albrechtsen et al., 2009) which allows SNPs
#' to be in LD and implicitly accounts for the LD through conditional emission probabilities where
#' the current genotype probability is conditioned on the genotype of a single previous SNP (haplotype frequencies).
#' }
#' \code{model=1} requires significantly less time to run than \code{model=2} due to the reduced number of SNPs
#' and simplified emission probabilities. It may be wise to use \code{model=1} when there are many SNPs and many
#' samples.
#' @param maf the smallest minor allele frequency allowed in the analysis. The
#' default value is 0.01.
#' @param sample.max.missing the maximum proportion of missing data allowed for
#' each sample. The default value is 0.1.
#' @param snp.max.missing the maximum proportion of missing data allowed for each
#' SNP. The default value is 0.1.
#' @param maximum.ld.r2 the maximum linkage disequilibrium R2 value allowed
#' between pairs of SNPs.The default value is 0.99.
#' @param chromosomes a numeric vector containing a subset of chromosomes to perform genotype filtering on. The
#' default is \code{chromosomes=NULL} which will format genotypes for all chromosomes in \code{ped.map}.
#' Autosomes are represented by numbers 1-22 and the X chromosome is denoted 23.
#' @param input.map.distance either "M" or "cM" denoting whether the genetic map distances in
#' \code{ped.map} are in Morgans (M) or centi-Morgans (cM). The default is Morgans.
#' @param reference.map.distance either "M" or "cM" denoting whether the genetic map distances in
#' \code{reference.ped.map} are in Morgans (M) or centi-Morgans (cM). The default is Morgans. HapMap reference data is in Morgans.
#' @return A named list of three objects:
#' \enumerate{
#' \item A pedigree containing the samples that remain after filtering. The pedigree is the first six columns
#' of the PED file and these columns are
#' headed \code{fid, iid, pid, mid, sex} and \code{aff}, respectively.
#' \item A data frame with the first five columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (Morgans, M) (type \code{"numeric"})
#' \item Base-pair position (type \code{"numeric"} or \code{"integer"})
#' \item Population allele frequency (type \code{"numeric"})
#' }
#' where each row describes a single marker. These columns are headed \code{chr, snp_id, pos_M, pos_bp} and \code{freq} respectively.
#' If \code{model=2} then the following columns are also included:
#' \enumerate{
#' \item Numeric ID of condition SNP (type \code{"numeric"} or \code{"integer"})
#' \item Haplotype probability: pba (type \code{"numeric"})
#' \item Haplotype probability: pbA (type \code{"numeric"})
#' \item Haplotype probability: pBa (type \code{"numeric"})
#' \item Haplotype probability: pBA (type \code{"numeric"})
#' \item Population allele frequency on the condition SNP (type \code{"numeric"})
#' }
#' with the headers \code{condition_snp, pba, pbA, pBa, pBA} and \code{freq_condition_snp}.
#' The remaining columns contain the genotype data for each sample, where a single column corresponds to a single sample. These columns are
#' labeled with merged family IDs and individual IDs separated by a slash symbol (/).
#' \item The model selected.
#' }
#' The list is named \code{pedigree, genotypes} and \code{model} respectively.
#' @export
#' @examples
#' # look at the simulated data
#' str(example_pedmap)
#'
#' # format and filter the example data using model 2 and reference data
#' my_genotypes <- getGenotypes(ped.map = example_pedmap,
#'                              reference.ped.map = example_reference_pedmap,
#'                              snp.ld = example_reference_ld,
#'                              model = 2,
#'                              maf = 0.01,
#'                              sample.max.missing = 0.1,
#'                              snp.max.missing = 0.1,
#'                              maximum.ld.r2 = 0.99,
#'                              chromosomes = NULL,
#'                              input.map.distance = "M",
#'                              reference.map.distance = "M")
getGenotypes <- function(ped.map, reference.ped.map = NULL, snp.ld = NULL, model = 1, maf = 0.01, sample.max.missing = 0.1,
                         snp.max.missing = 0.1, maximum.ld.r2 = 0.99, chromosomes = NULL, input.map.distance = "M",
                         reference.map.distance = "M"){

  # check the input parameters

  # check input PED and MAP files
  stopifnot(is.list(ped.map) | length(ped.map) == 2)
  input.ped <- ped.map[[1]]
  input.map <- ped.map[[2]]

  # check the PED and MAP files have the same number of SNPs
  if (ncol(input.ped) != (2*nrow(input.map)+6))
    stop ("PED and MAP files are not in the correct format")

  # check the MAP file has 4 coloumns
  if (ncol(input.map) != 4)
    stop ("MAP file has incorrect format")
  colnames(input.map) <- c("chr", "snp_id", "pos_M","pos_bp")
  if (!is.integer(input.map[,"chr"]))
    stop ("chromosomes in MAP file must be numeric")

  # check reference data
  if (!is.null(reference.ped.map)) {
    stopifnot(is.list(reference.ped.map) & length(reference.ped.map) == 2)
    reference.ped <- reference.ped.map[[1]]
    reference.map <- reference.ped.map[[2]]

    # check the PED and MAP files have the same number of SNPs
    if (ncol(reference.ped) != (2*nrow(reference.map)+6))
      stop ("reference PED and MAP files are not in the correct format")

    # check the MAP file has 4 coloumns
    if (ncol(reference.map) != 4)
      stop ("reference MAP file has incorrect format")
    colnames(reference.map) <- c("chr", "snp_id", "pos_M","pos_bp")
  }

  # check snp.ld
  if (!is.null(snp.ld)) {
    stopifnot(ncol(snp.ld) == 7 & colnames(snp.ld) == c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"))
    if (model != 2)  {
      warning(paste0("'snp.ld' will be ignored since 'model=",model,"'"))
      snp.ld <- NULL
    }
  } else {
    if (model == 2)
      stop("'model=2' requires 'snp.ld'")
  }

  # check numeric input parameters
  stopifnot(is.numeric(maf))
  stopifnot(is.numeric(sample.max.missing))
  stopifnot(is.numeric(snp.max.missing))
  stopifnot(is.numeric(model))
  if (model != 1 & model != 2) {
    cat(paste("'model=",model,"' is not valid. Setting 'model=1'",sep=""))
    model <- 1
  }

  # check chromosomes
  if (!is.null(chromosomes)) {
    stopifnot(is.numeric(chromosomes))
    if(!all(chromosomes %in% input.map[,"chr"]))
      stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% input.map[,"chr"])]," not in 'ped.map'\n")))
  } else
    chromosomes <- as.numeric(unique(as.character(input.map[,"chr"])))
  chromosomes <- chromosomes[order(chromosomes)]

  # check input map distance
  if (input.map.distance != "M" & input.map.distance != "cM")
    stop (paste0("'input.map.distance=",input.map.distance,"' is not valid"))
  if (input.map.distance == "cM") {
    input.map[,"pos_M"] <- input.map[,"pos_M"]/100
  }

  # check reference map distance
  if (reference.map.distance != "M" & reference.map.distance != "cM")
    stop (paste0("'reference.map.distance=",reference.map.distance,"' is not valid"))
  if (reference.map.distance == "cM") {
    reference.map[,"pos_M"] <- reference.map[,"pos_M"]/100
  }

  # begin data filtering
  cat(paste("Begin filtering of ",nrow(input.ped)," samples and ",nrow(input.map)," SNPs...\n",sep=""))

  # create new sample IDs from PED FIDs and IIDs
  sample.names <- paste(input.ped[,1], input.ped[,2], sep="/")


  # merge input data with reference data & subset by chromosomes
  if (!is.null(reference.ped.map)) {
    input.map.v1      <- cbind(1:nrow(input.map), input.map)
    reference.map.v1  <- cbind(1:nrow(reference.map), reference.map)
    input.map.v1      <- merge(input.map.v1, reference.map.v1, by.x="snp_id", by.y="snp_id")
    input.map.v1      <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
    if (!is.null(chromosomes))
      input.map.v1    <- input.map.v1[input.map.v1[,"chr.x"] %in% chromosomes,]
    input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
    input.ped.columns <- input.ped.columns[order(input.ped.columns)]
    input.ped.v1      <- input.ped[,input.ped.columns]
    reference.ped.columns  <- c(1:6, 2*input.map.v1[,"1:nrow(reference.map)"] + 5, 2*input.map.v1[,"1:nrow(reference.map)"] + 6)
    reference.ped.columns  <- reference.ped.columns[order(reference.ped.columns)]
    reference.ped.v1       <- reference.ped[,reference.ped.columns]
    input.map.v2           <- input.map.v1[,c("chr.y", "snp_id", "pos_M.y", "pos_bp.y")]
    colnames(input.map.v2) <- c("chr", "snp_id", "pos_M","pos_bp")
  } else {
    input.map.v1      <- cbind(1:nrow(input.map), input.map)
    input.map.v1      <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
    if (!is.null(chromosomes))
      input.map.v1    <- input.map.v1[input.map.v1[,"chr"] %in% chromosomes,]
    input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
    input.ped.columns <- input.ped.columns[order(input.ped.columns)]
    input.ped.v1      <- input.ped[,input.ped.columns]
    input.map.v2      <- input.map.v1[,c("chr", "snp_id", "pos_M", "pos_bp")]
  }
  if (nrow(input.map.v2) == 0)
    stop("0 SNPs remain after merging with reference dataset")
  cat(paste(nrow(input.map.v2)," SNPs remain after merging with reference dataset...\n",sep=""))


  # call genotypes
  input.matrix        <- as.matrix(input.ped.v1[,7:ncol(input.ped.v1)])
  input.genders       <- input.ped.v1[,5]
  input.chromosomes   <- input.map.v2[,1]
  input.genotypes.v0  <- cbind(input.map.v2, haplotypeToGenotype(input.matrix, input.chromosomes, input.genders))
  if (!is.null(reference.ped.map)) {
    reference.matrix       <- as.matrix(reference.ped.v1[,7:ncol(reference.ped.v1)])
    reference.genders      <- reference.ped.v1[,5]
    reference.genotypes.v0 <- cbind(input.map.v2, haplotypeToGenotype(reference.matrix, input.chromosomes, reference.genders))
  }


  # calculate allele frequencies form reference data
  if (is.null(reference.ped.map)) {
    pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(input.genotypes.v0[,5:ncol(input.genotypes.v0)]), as.numeric(input.genotypes.v0[,"chr"]), input.ped.v1[,5])
    input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,5:ncol(input.genotypes.v0)])
  } else {
    pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(reference.genotypes.v0[,5:ncol(reference.genotypes.v0)]), as.numeric(reference.genotypes.v0[,"chr"]), reference.ped.v1[,5])
    input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,c(5:ncol(input.genotypes.v0))])
    if (model == 2){
      reference.genotypes.v1 <- cbind(reference.genotypes.v0[,c(1:4)],pop.allele.freq,reference.genotypes.v0[,c(5:ncol(reference.genotypes.v0))])
      colnames(reference.genotypes.v1)[1:5] <- c("chr", "snp_id", "pos_M","pos_bp", "freq")
    }
  }
  colnames(input.genotypes.v1) <- c("chr", "snp_id", "pos_M","pos_bp", "freq", sample.names)
  #cat(paste("Begin filtering of ",length(sample.names)," samples and ",nrow(input.genotypes.v1)," SNPs...\n",sep=""))


  # remove SNPs with low population MAF
  input.genotypes.v2 <- subset(input.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf)
  if (nrow(input.genotypes.v2) == 0)
    stop("0 SNPs remain after MAF removal")
  if (!is.null(reference.ped.map) & model == 2)
    reference.genotypes.v2 <- subset(reference.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf)
  cat(paste(nrow(input.genotypes.v2)," SNPs remain after MAF removal...\n",sep=""))


  # remove snps with high missingness
  snp.missingness <- calculateMissingness(as.matrix(t(input.genotypes.v2[,6:ncol(input.genotypes.v2)])))
  if (!is.null(reference.ped.map) & model == 2) {
    snp.missingness.ref    <- calculateMissingness(as.matrix(t(reference.genotypes.v2[,6:ncol(reference.genotypes.v2)])))
    input.genotypes.v3     <- input.genotypes.v2[snp.missingness <= snp.max.missing & snp.missingness.ref <= snp.max.missing,]
    reference.genotypes.v3 <- reference.genotypes.v2[snp.missingness <= snp.max.missing & snp.missingness.ref <= snp.max.missing,]
  } else {
    input.genotypes.v3 <- input.genotypes.v2[snp.missingness <= snp.max.missing,]
  }
  if (nrow(input.genotypes.v3) == 0)
    stop("0 SNPs remain after missingness removal")
  cat(paste(nrow(input.genotypes.v3)," SNPs remain after missingness removal...\n",sep=""))


  # remove SNPs in high LD
  if (model == 2){
    highLD <- c(snp.ld[snp.ld[,"R2"] > maximum.ld.r2, "SNP_A"], snp.ld[snp.ld[,"R2"] > maximum.ld.r2, "SNP_B"])
    if (length(highLD) > 0) {
      highLD <- unique(highLD)
      input.genotypes.v4 <- input.genotypes.v3[!(input.genotypes.v3[,"snp_id"] %in% highLD),]
      #input.genotypes.v4 <- subset(input.genotypes.v3, !(snp_id %in% highLD))
      if (!is.null(reference.ped.map) & model == 2)
        reference.genotypes.v4 <- reference.genotypes.v3[!(reference.genotypes.v3[,"snp_id"] %in% highLD),]
        #reference.genotypes.v4 <- subset(reference.genotypes.v3, !(snp_id %in% highLD))
    } else {
      input.genotypes.v4 <- input.genotypes.v3
      if (!is.null(reference.ped.map) & model == 2)
        reference.genotypes.v4 <- reference.genotypes.v3
    }
    if (nrow(input.genotypes.v4) == 0)
      stop("0 SNPs remain after LD removal")
    cat(paste(nrow(input.genotypes.v4)," SNPs remain after LD removal...\n",sep=""))
  } else {
    input.genotypes.v4 <- input.genotypes.v3
    if (!is.null(reference.ped.map) & model == 2)
      reference.genotypes.v4 <- reference.genotypes.v3
  }


  # remove samples with high missingness
  sample.missingness <- round(calculateMissingness(as.matrix(input.genotypes.v4[,6:ncol(input.genotypes.v4)])),digits=3)
  if (length(sample.names[sample.missingness > sample.max.missing]) > 0) {
    cat(paste("*********\n WARNING: Removing ",sample.names[sample.missingness > sample.max.missing]," because genotype missingness > ",sample.max.missing*100,"%\n*********\n",sep=""))
    sample.keep        <- input.ped.v1[sample.missingness <= sample.max.missing,1:6]
    input.genotypes.v5 <- input.genotypes.v4[,c(1:5, which(sample.missingness <= sample.max.missing) + 5)]
    if(nrow(sample.keep) < 1) stop(paste("All samples removed with missingness > ",sample.max.missing*100,"%. No samples remaining.",sep=""))
  } else {
    sample.keep        <- input.ped.v1[,1:6]
    input.genotypes.v5 <- input.genotypes.v4
  }
  colnames(sample.keep) <- c("fid", "iid", "pid", "mid", "sex", "aff")
  if ((ncol(input.genotypes.v5)-5) == 0) {
    stop("0 samples remain after missingness removal")
  }
  cat(paste(ncol(input.genotypes.v5)-5," samples remain after missingness removal...\n",sep=""))


  # find SNPs to condition on and haplotype frequencies
  if (model == 2){
    if (is.null(reference.ped.map)) {
      genders <- sample.keep[,5]
      condition.snps <- getConditionSNP(input.genotypes.v5, snp.ld, maximum.ld.r2)
      condition.matrix <- getHaplotypeFreq(input.genotypes.v5, condition.snps, genders)
    } else {
      genders <- reference.ped.v1[,5]
      condition.snps <- getConditionSNP(reference.genotypes.v4, snp.ld, maximum.ld.r2)
      condition.matrix <- getHaplotypeFreq(reference.genotypes.v4, condition.snps, genders)
    }
    input.genotypes.v6 <- cbind(condition.matrix, input.genotypes.v5[,6:ncol(input.genotypes.v5)])
  } else {
    input.genotypes.v6 <- input.genotypes.v5
  }


  return.genotypes <- list(sample.keep, input.genotypes.v6, model)
  names(return.genotypes) <- c("pedigree", "genotypes", "model")
  return(return.genotypes)
}

