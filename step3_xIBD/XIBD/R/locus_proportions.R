#' XIBD Proportion of Pairs IBD
#'
#' Calculate the proportion of pairs IBD at each SNP.
#'
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param locus.matrix a data frame containing the binary IBD information for each SNP and each pair.
#' See \code{value} description in \code{\link{getLocusMatrix}} for more details.
#' @param groups a data frame with three columns of information:
#' \enumerate{
#' \item Family ID
#' \item Sample ID
#' \item Group ID
#' }
#' Group ID, for example, can be affection status.
#' If \code{groups} is specified, each sample in the pedigree should belong to a group,
#' and locus proportions will be calculated for each group.
#' The default is \code{groups=NULL}.
#' @return A data frame the following 5 columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' \item Proportion of pairs IBD (type \code{"integer"})
#' }
#' where each row describes a unique SNP. The data frame is headed
#' \code{chr, snp_id, pos_M, pos_bp} and \code{prop} respectively. If \code{groups} is
#' specified then there will be one column of proportions of each combination of groups with
#' these column header by group IDs and the number of pairs in the group.
#' @export
#' @examples
#' # generate a binary IBD matrix
#' my_locus_matrix <- getLocusMatrix(ped.genotypes = example_genotypes,
#'                                   ibd.segments = example_ibd)
#'
#' # calculate the proportion of pairs IBD at each SNP
#' my_locus_prop <- getLocusProportion(ped.genotypes = example_genotypes,
#'                                     locus.matrix = my_locus_matrix,
#'                                     groups = NULL)
getLocusProportion <- function(ped.genotypes, locus.matrix, groups = NULL){

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 2)
  pedigree <- ped.genotypes[[1]]

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check locus matrix input
  if (ncol(locus.matrix) < 5)
    stop ("locus.matrix has incorrect format")
  colnames(locus.matrix)[1:4] <- c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp")

  # check groups
  if (!is.null(groups)) {
    stopifnot(is.data.frame(groups))
    stopifnot(ncol(groups) > 2)
    if (ncol(groups) > 3){
      cat("using first 3 columns of groups")
      groups <- groups[,1:3]
    }
    colnames(groups)[1:2] <- c("fid","iid")

    # check isolates belong to a group
    group.names <- paste(groups[,"fid"],groups[,"iid"],sep="/")
    isolate.names <- paste(pedigree[,"fid"],pedigree[,"iid"],sep="/")
    if (!all(isolate.names %in% group.names))
      stop("'groups' is missing information for some isoaltes")

    # assign number ID to each isoalte
    pedigree.0 <- data.frame(num.id=1:nrow(pedigree), pedigree)

    # merge pedigree with proups by IDs
    pedigree.group <- merge(pedigree.0, groups, by=c("fid", "iid"))

    # reorder merged pedigree by numberic IDs
    pedigree.group <- pedigree.group[order(pedigree.group[,"num.id"]),]

    # get isolates group pairs - dataframe with 2 coloumns
    group.pairs <- groupPairs(as.character(pedigree.group[,8]))

    # reorder pairs groups
    groups.unique <- as.character(unique(pedigree.group[,8]))
    groups.unique.pairs <- groupPairs(groups.unique)
    if (nrow(groups.unique.pairs) > 0) {
      for (i in 1:nrow(groups.unique.pairs)){
        change.pair <- which(group.pairs[,1] == groups.unique.pairs[i,2] & group.pairs[,2] == groups.unique.pairs[i,1])
        group.pairs[change.pair,1] <- groups.unique.pairs[i,1]
        group.pairs[change.pair,2] <- groups.unique.pairs[i,2]
      }
    }

    # get number of pairwise analyses for each group
    number.pairs.groups <- NULL
    for (i in unique(paste(group.pairs[,1],group.pairs[,2],sep="/"))) {
      id_1 <- unlist(strsplit(i,"/"))[1]
      id_2 <- unlist(strsplit(i,"/"))[2]
      number.pairs.groups[group.pairs[,1] == id_1 & group.pairs[,2] == id_2] <- paste("npairs=",
                                                                                      dim(group.pairs[group.pairs[,1] == id_1 & group.pairs[,2] == id_2,])[1],
                                                                                      sep="")
    }

    # combine group pairs into a single vector
    group.pairs.1 <- paste(group.pairs[,1],group.pairs[,2],number.pairs.groups,sep="\n")
  }

  # calculate proportion IBD at each SNP
  locus.pairs <- locus.matrix[,5:ncol(locus.matrix)]
  if (!is.null(groups)) {
    locus.prop <- NULL
    for (g in unique(group.pairs.1)) {
      locus.pairs.g <- locus.pairs[,group.pairs.1 == g]
      locus.prop <- cbind(locus.prop, rowMeans(locus.pairs.g))
    }
    colnames(locus.prop) <- unique(group.pairs.1)
  } else {
    locus.prop <- rowMeans(locus.pairs)
  }

  return.locus.prop <- cbind(locus.matrix[,1:4],locus.prop)
  if(is.null(groups))
    colnames(return.locus.prop)[5] <- "prop"

  return(data.frame(return.locus.prop))
}


