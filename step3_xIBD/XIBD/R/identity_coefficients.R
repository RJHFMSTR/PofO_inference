#' to prevent notes
globalVariables("pair.i")
#' Identity Coefficients
#'
#' Determines the identity coefficients based on a pedigree.
#' Autosome coefficients are calculated using the R package \code{identity} while
#' the X chromosome coefficients are implemented by \code{XIBD}.
#'
#'
#' @param pedigree a data frame containing a pedigree with 6 columns
#' \enumerate{
#' \item Family ID (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item Individual ID (\strong{\code{"numeric"}} or \strong{\code{"integer"}})
#' \item Paternal ID (\strong{\code{"numeric"}} or \strong{\code{"integer"}})
#' \item Maternal ID (\strong{\code{"numeric"}} or \strong{\code{"integer"}})
#' \item Gender (1 = male, 2 = female)
#' \item Phenotype (1 = unaffected, 2 = affected, 0 = unknown)
#' }
#' The family ID is not actually used in the analysis but is required for completeness of the pedigree;
#' thus \strong{all individual IDs must be unique and numeric.} \strong{The individuals in the pedigree
#' should also be numbered in a way such that every parent precedes his or her children.}
#' @param number.cores the number of cores used for parallel execution.
#' @return A data frame with columns:
#' \enumerate{
#' \item Family 1 ID (type \code{"character"})
#' \item Individual 1 ID (type \code{"character"})
#' \item Family 2 ID (type \code{"character"})
#' \item Individual 2 ID (type \code{"character"})
#' \item The number of meiosis (type \code{"numeric"} or \code{"integer"})
#' \item Autosome identity coefficient for sharing 0 alleles IBD (type \code{"numeric"})
#' \item Autosome identity coefficient for sharing 1 allele IBD (type \code{"numeric"})
#' \item Autosome identity coefficient for sharing 2 alleles IBD (type \code{"numeric"})
#' \item X chromsome identity coefficient for sharing 0 alleles IBD (type \code{"numeric"})
#' \item X chromsome identity coefficient for sharing 1 allele IBD (type \code{"numeric"})
#' \item X chromsome identity coefficient for sharing 2 alleles IBD (type \code{"numeric"})
#' }
#' headed \code{fid1, iid1, fid2, iid2, a_i0, a_i1, a_i2, x_i0, x_i1} and \code{x_i2}, respectively.
#' Each row describes identity coefficient for a unique pair of samples.
#' @importFrom foreach "%dopar%"
#' @importFrom identity identity.coefs
#' @export
#' @examples
#' \dontrun{
#' # generate a pedigree
#' my_pedigree <- data.frame(fid = rep(1,16),
#'                           iid = 1:16,
#'                           pid = c(0,0,0,1,1,0,0,4,6,0,7,0,0,9,11,13),
#'                           mid = c(0,0,0,2,2,0,0,3,5,0,8,0,0,10,12,14),
#'                           sex = c(1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1),
#'                           aff = rep(1,16))
#'
#' # calculate identity coefficients
#' identity_coef <- getIdentityCoef(pedigree = my_pedigree,
#'                                  number.cores = 1)
#'
#' head(identity_coef)
#' }
getIdentityCoef <- function(pedigree, number.cores) {
  # check pedigree
  stopifnot(ncol(pedigree) == 6)
  stopifnot(nrow(pedigree) > 2)
  if (!is.numeric(pedigree[,"iid"]))
    stop("'pedigree[,iid]' is not numeric")
  if (!is.numeric(pedigree[,"pid"]))
    stop("'pedigree[,pid]' is not numeric")
  if (!is.numeric(pedigree[,"mid"]))
    stop("'pedigree[,mid]' is not numeric")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  stopifnot(is.numeric(number.cores))

  # get sample pairs
  sample_pairs <- isolatePairs(pedigree[,"fid"],pedigree[,"iid"])

  # turn pedigree into a matrix for
  new_pedigree <- as.matrix(pedigree[,c("iid", "pid", "mid", "sex")])


  # set number of cores for R package parallel
  doParallel::registerDoParallel(cores=number.cores)

  # calculating coefficients
  ibd_coefs <- foreach::foreach(pair.i=1:nrow(sample_pairs), .combine='rbind') %dopar% {
    fid.1 <- as.character(sample_pairs[pair.i,1])
    ind.1 <- as.numeric(sample_pairs[pair.i,2])
    fid.2 <- as.character(sample_pairs[pair.i,3])
    ind.2 <- as.numeric(sample_pairs[pair.i,4])
    a.coef <- matrix(identity::identity.coefs(c(ind.1, ind.2), new_pedigree)[2,],nrow=1,ncol=11)
    ibd.coefs <- cbind(fid.1, ind.1, fid.2, ind.2,
                       ibd_coefs_A(a.coef),
                       ibd_coefs_X(new_pedigree, ind.1, ind.2))
    ibd.coefs
  }
  colnames(ibd_coefs) <- c("fid1","iid1","fid2","iid2","a_i0","a_i1","a_i2","x_i0","x_i1","x_i2")

  # change structure of data
  ibd_coefs <- data.frame(ibd_coefs)
  for(i in 1:4){
    ibd_coefs[,i] <- as.character(ibd_coefs[,i])
  }
  for(i in 5:10){
    ibd_coefs[,i] <- as.numeric(as.character(ibd_coefs[,i]))
  }

  return(ibd_coefs)
}






