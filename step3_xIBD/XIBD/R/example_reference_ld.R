#' Reference Linkage Disequilibrium Data for Example
#'
#'
#' A subset of the HapMap TSI linkage disequilibrium (LD) data for SNPs
#' in the \code{example_pedmap.rda} file.
#'
#' @format Data frame with columns
#' \describe{
#' \item{CHR_A}{Chromosomes of SNP A.}
#' \item{BP_A}{Base-pair positions of SNP A.}
#' \item{SNP_A}{SNP identifiers of SNP A.}
#' \item{CHR_B}{Chromosomes of SNP B.}
#' \item{BP_B}{Base-pair positions of SNP B.}
#' \item{SNP_B}{SNP identifiers of SNP B.}
#' \item{R2}{Squared correlation between SNPs A and B.}
#' }
#' See \code{\link{getGenotypes}} for details.
"example_reference_ld"
