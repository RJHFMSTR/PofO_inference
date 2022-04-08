#' IBD Segments For The Simulated Data
#'
#'
#' The IBD segments inferred using XIBD with the
#' parameter settings as in the Vignette.
#'
#' @format A data frame with
#'  \describe{
#' \item{fid1}{Family 1 ID}
#' \item{iid1}{Individual 1 ID}
#' \item{fid2}{Family 2 ID}
#' \item{iid2}{Individual 2 ID}
#' \item{chr}{Chromosome}
#' \item{start_snp}{IBD segment start SNP ID}
#' \item{end_snp}{IBD segment end SNP ID}
#' \item{start_position_bp}{IBD segment start SNP base-pair position}
#' \item{end_position_bp}{IBD segment end SNP base-pair position}
#' \item{start_position_M}{IBD segment start SNP morgan position}
#' \item{end_position_M}{IBD segment end SNP morgan position}
#' \item{number_snps}{Number of SNPs in IBD segment}
#' \item{length_bp}{Length of IBD segment in base-pairs}
#' \item{length_M}{Length of IBD segment in morgans}
#' \item{ibd_status}{The number of alleles shared IBD (either 1 or 2)}
#' }
"example_ibd"
