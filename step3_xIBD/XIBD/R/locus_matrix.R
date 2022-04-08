#' XIBD Binary IBD Matrix
#'
#' Returns a binary matrix of IBD (1) and non-IBD (0) for each SNP and pair combination.
#'
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param ibd.segments a named list containing the \code{ibd_segments} inferred from pairs of samples.
#' See \code{value} description in \code{\link{getIBDsegments}} for more details.
#' @return A data frame the first four columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' }
#' where each row describes a unique SNP. These columns are headed \code{chr, snp_id, pos_M} and \code{pos_bp} respectively.
#' Columns 5 onwards contain the binary IBD information for each sample pair, where a single column corresponds to a single pair.
#' These columns are labeled with merged family IDs and individual IDs separated by a slash symbol (/). For example famA/ind1/famA/ind2.
#' @export
#' @examples
#' # generate a binary IBD matrix
#' my_locus_matrix <- getLocusMatrix(ped.genotypes = example_genotypes,
#'                                   ibd.segments = example_ibd)
getLocusMatrix <- function(ped.genotypes, ibd.segments){

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  if(!all(c("pedigree", "genotypes", "model") %in% names(ped.genotypes)))
    stop("the objects 'pedigree', 'genotypes' and 'model' are not in 'ped.genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) < 8 & nrow(genotypes) <= 1)
    stop ("ped.genotypes has incorrect format")
  colnames(genotypes)[1:5] <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp", "FREQ")

  # check ibd.segments file is a dataframe with correct fields
  stopifnot(is.list(ibd.segments))
  stopifnot("ibd_segments" %in% names(ibd.segments))
  stopifnot(is.data.frame(ibd.segments[["ibd_segments"]]))
  ibd.segments <- ibd.segments[["ibd_segments"]]
  if (ncol(ibd.segments) != 15)
    stop ("ibd.segments has incorrect format")
  colnames(ibd.segments) <- c("fid1", "ind1", "fid2", "ind2", "chr", "start.snp", "end.snp", "start.position.bp",
                              "end.position.bp", "start.position.M", "end.position.M", "number.snps", "length.bp",
                              "length.M", "ibd.status")

  # create a data frame of pairs; one unique pair per row; columns are fid1, iid1, fid2, iid2 and
  # a unique numeric pair identifier
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  isolate.pairs <- cbind(isolate.pairs, 1:nrow(isolate.pairs))
  colnames(isolate.pairs) <- c("fid1","ind1","fid2","ind2","id")

  # merge inferred IBD pairs with their numeric ID
  isolate.id <- merge(isolate.pairs, ibd.segments)

  # generate the binary IBD matrix
  ibd.binary.matrix <- IBDMatrix(as.character(genotypes[,"CHROMOSOME"]), genotypes[,"POSITION.bp"], nrow(isolate.pairs),
                                 as.numeric(as.character(isolate.id[,"id"])), isolate.id[,"chr"], isolate.id[,"start.position.bp"],
                                 isolate.id[,"end.position.bp"])
  ibd.locus.matrix <- data.frame(genotypes[,c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp")],ibd.binary.matrix)
  colnames(ibd.locus.matrix) <- c("chr","snp_id","pos_M","pos_bp",
                                  paste(isolate.pairs[,"fid1"],isolate.pairs[,"ind1"],
                                        isolate.pairs[,"fid2"],isolate.pairs[,"ind2"], sep="/"))

  return(data.frame(ibd.locus.matrix))
}


