#' Filtered Genotypes For The Simulated Data
#'
#'
#' Processed raw genotype data with the
#' parameter settings as in the Vignette
#' and following topbot alterations.
#'
#' @format A list of two objects named \code{pedigree} and \code{genotypes}:
#' \enumerate{
#' \item \code{pedigree} is a data frame with the following information:
#' \enumerate{
#' \item{Family ID}
#' \item{Isolate ID}
#' \item{Paternal ID. This is not used by XIBD and is set to zero.}
#' \item{Maternal ID. This is not used by XIBD and is set to zero.}
#' \item{Gender - 1 for male and 2 for female}
#' \item{Affection status of individuals. This is set to 2 and is ignored by XIBD IBD inference.}
#' }
#'
#' \item \code{genotypes} is a data frame with the first 11 columns:
#' \enumerate{
#' \item Chromosome
#' \item SNP identifier
#' \item Genetic map distance
#' \item Base-pair position
#' \item Population allele frequency
#' \item Condition SNP number
#' \item Conditional probability Pr(ba)
#' \item Conditional probability Pr(bA)
#' \item Conditional probability Pr(Ba)
#' \item Conditional probability Pr(BA)
#' \item Allele frequency of conditioned SNP
#' }
#' where each row describes a single SNP. Columns 12 onwards contain the genotype data for each individual, where a single column corresponds to a single individual.
#' These columns are labeled with merged family IDs and individual IDs separated by a slash symbol (/).
#' }
"example_genotypes"
