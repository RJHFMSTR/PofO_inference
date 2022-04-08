#' XIBD Calculate Population Allele Frequencies
#'
#' Calculate population allele frequencies from genotypes.
#'
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' @return a data frame with information:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (Morgans, M) (type \code{"numeric"})
#' \item Base-pair position (type \code{"numeric"} or \code{"integer"})
#' \item Population allele frequency (type \code{"numeric"})
#' }
#' where each row defines a unique SNP.
#' Allele frequencies of -1 indicate missing data for all samples for that SNP.
#' @export
calculateAlleleFreq <- function(ped.genotypes) {

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  stopifnot(c("pedigree","genotypes","model") %in% names(ped.genotypes))
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  model     <- ped.genotypes[["model"]]
  stopifnot(is.numeric(model))
  stopifnot(model == 1 | model == 2)

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("'ped.genotypes' has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (model == 1){
    if(ncol(genotypes) < 8 & nrow(genotypes) <= 1)
      stop ("'ped.genotypes' has incorrect format")
    if(!all(colnames(genotypes)[1:5] %in% c("chr", "snp_id", "pos_M","pos_bp", "freq")))
      stop ("'ped.genotypes' has incorrect format")
    samples <- colnames(genotypes)[!(colnames(genotypes) %in% c("chr", "snp_id", "pos_M","pos_bp", "freq"))]
  }
  if (model == 2) {
    if(ncol(genotypes) < 14 & nrow(genotypes) <= 1)
      stop ("'ped.genotypes' has incorrect format")
    if(!all(colnames(genotypes)[1:11] %in% c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")))
      stop ("'ped.genotypes' has incorrect format")
    samples <- colnames(genotypes)[!(colnames(genotypes) %in% c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp"))]
  }

  my.chr <- genotypes[,"chr"]
  my.sex <- pedigree[,"sex"]
  my.geno <- as.matrix(genotypes[,paste(pedigree[,"fid"], pedigree[,"iid"], sep="/")],
                       nrow=length(my.chr), ncol=length(my.sex))
  my.freq <- calculatePopAlleleFreq(my.geno, my.chr, my.sex)
  return.freq <- data.frame(genotypes[,c("chr","snp_id","pos_M","pos_bp")],freq=my.freq)

  return(return.freq)
}


