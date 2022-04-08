# Internal Function
#
# Merge Returned Lists from Parallele
#
# \code{merge_lists()} is a function used to merge summary IBD results for multiple pairs when running the IBD analysis on
# multiple cores
#
# @param A List with n objects for one pair.
# @param B List with n objects for another pair. The dimension of each object in A and B should equal.
# @return A list with 2 objects containing merged lists from \code{A} and \code{B}    above.
mergeLists1 <- function(A, B){
  x <- list()
  for(i in 1:max(length(A),length(B))) {
    x[[i]] <- mapply(rbind, A[i], B[i], SIMPLIFY=FALSE)[[1]]
  }
  return(x)
}


mergeLists2 <- function(A, B){
  x <- list()
  x[[1]] <- mapply(rbind, A[1], B[1], SIMPLIFY=FALSE)[[1]]
  x[[2]] <- mapply(cbind, A[2], B[2], SIMPLIFY=FALSE)[[1]]
  return(x)
}


#' @useDynLib XIBD
#' @importFrom Rcpp sourceCpp
NULL
