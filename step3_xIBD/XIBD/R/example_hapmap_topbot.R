#' Example Topbot Data
#'
#' A topbot annotation file containing Illumina TOP/BOT designation for HapMap SNPs. This is a subset of HapMap SNPs.
#'
#' @format Data frame containing
#' \describe{
#' \item{chr}{Chromosome (\code{"integer"}).}
#' \item{snp_id}{SNP identifier (type \code{"character"}).}
#' \item{pos_bp}{Base-pair position (type \code{"numeric"} or \code{"integer"}).}
#' \item{pos_M}{Genetic map distance (centi morgans, cM, or morgans, M) (type \code{"numeric"}).}
#' \item{TOPBOT}{Illuminas TOP or BOT designation of the SNP (type \code{"character"}).}
#' }
"example_hapmap_topbot"
