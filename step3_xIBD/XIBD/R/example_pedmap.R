#' Simulated IBD Data for Example
#'
#' Simulated haplotype data in PLINK PED/MAP format for 28,808 SNPs
#' from chromosomes 1-23.
#' The data is for 10 samples from a 5 generation family with the most
#' distant pair being third cousins and the shortest distant pair being
#' parent-offspring.
#' The data was simulated for SNPs from an Illumina SNPchip under Illuminas
#' TOPBOT allele naming convention (see \code{switchBOTgenotypes}).
#'
#' @format List with two objects
#' \describe{
#' \item{PED}{PLINK PED data frame. See \code{\link{getGenotypes}} for details.}
#' \item{MAP}{PLINK MAP data frame. See \code{\link{getGenotypes}} for details.}
#' }
"example_pedmap"
