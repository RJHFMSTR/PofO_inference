#' IBD Parameter Estimates For Simulated Data
#'
#'
#' The parameter estimates inferred using XIBD with the
#' parameter settings as in the Vignette.
#'
#' @format A list of two objects named \code{autosome_parameters} and \code{X_chromosome_parameters},
#' where each object contains the following fields.
#'  \describe{
#' \item{fid1}{Family 1 ID}
#' \item{iid1}{Individual 1 ID}
#' \item{fid2}{Family 2 ID}
#' \item{iid2}{Individual 2 ID}
#' \item{m}{The estimated number of meiosis separating each pair of individuals}
#' \item{ibd0}{The estimated proportion of genome with 0 alleles IBD}
#' \item{ibd1}{The estimated proportion of genome with 1 allele IBD}
#' \item{ibd2}{The estimated proportion of genome with 2 alleles IBD}
#' }
"example_parameters"
