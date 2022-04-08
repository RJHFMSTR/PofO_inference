# Internal Function
#
# IBD Summary Table
#
# \code{IBDTable()} produces summaries of detected IBD segments for a single pair of isolates. These summaries include the genetic
# map start and end
# of IBD segments in bp, cM and SNP identifiers; and lengths of IBD segments in bp, cM and SNPs.
# @param ibd.results A data frame containing family ID and isolate ID for isolate 1, family ID and isolate ID for isolate 2,
# numeric SNP identifiers, chromosome identifiers, genetic map positions of SNPs in Morgans (M) and base-pairs (bp) and the Viterbi results respectively,
# for all SNPs.
# @return A data frame containing a summary of all IBD segments inferred for this pair of isolates. The data frame contains
# the following columns:
# \enumerate{
# \item Family 1 ID
# \item Isolate 1 ID
# \item Family 2 ID
# \item Isolate 2 ID
# \item Chromosome
# \item SNP identifier
# \item Start SNP
# \item End SNP
# \item Start position bp
# \item End position bp
# \item Start position M
# \item End position M
# \item Number of SNPs
# \item Length bp
# \item Length M
# \item IBD status (1 = 1 allele shared IBD, 2 = 2 alleles shared IBD)
# }
IBDTable <- function(ibd.results){

  # turn the IBD results into a matrix (seems to be the only way this code works)
  ibd.colnames <- c("fid1","iid1","fid2","iid2","markerNo","chr","marker","pos.m","pos.bp","viterbi")
  ibd.matrix <- as.matrix(ibd.results, ncol=10)
  colnames(ibd.matrix) <- ibd.colnames

  # subset IBD matrix by either IBD=1 or IBD=2 segments
  ibd.table <- NULL
  for (IBD in 1:2) {

    ibd.positive.1 <- matrix(ibd.matrix[ibd.matrix[,"viterbi"] == IBD,], ncol=10)
    colnames(ibd.positive.1) <- ibd.colnames

    if (dim(ibd.positive.1)[1] != 0) {

      # subset results by chromosome
      for (chrom in unique(ibd.positive.1[,"chr"])) {

        ibd.positive.2 <- matrix(ibd.positive.1[ibd.positive.1[,"chr"] == chrom,], ncol=10)
        colnames(ibd.positive.2) <- ibd.colnames

        # number each distinct IBD segment
        ibd.positive.3 <- cbind(IBDLabel(as.numeric(ibd.positive.2[,"markerNo"]), dim(ibd.positive.2)[1]), ibd.positive.2)
        colnames(ibd.positive.3)[1] <- "ibd.label"

        # for each unique IBD segment, set summary info
        for (tract in unique(ibd.positive.3[,"ibd.label"])) {

          ibd.positive.4             <- matrix(ibd.positive.3[ibd.positive.3[,"ibd.label"] == tract,], ncol=11)
          colnames(ibd.positive.4)   <- colnames(ibd.positive.3)
          fid1              <- unique(ibd.positive.4[,"fid1"])
          fid2              <- unique(ibd.positive.4[,"fid2"])
          iid1              <- unique(ibd.positive.4[,"iid1"])
          iid2              <- unique(ibd.positive.4[,"iid2"])
          chr               <- unique(ibd.positive.4[,"chr"])
          number.snps       <- dim(ibd.positive.4)[1]
          start.position.bp <- as.numeric(ibd.positive.4[1,"pos.bp"])
          end.position.bp   <- as.numeric(ibd.positive.4[number.snps,"pos.bp"])
          start.position.M  <- as.numeric(ibd.positive.4[1,"pos.m"])
          end.position.M    <- as.numeric(ibd.positive.4[number.snps,"pos.m"])
          length.bp         <- end.position.bp - start.position.bp
          length.M          <- as.numeric(ibd.positive.4[number.snps,"pos.m"]) - as.numeric(ibd.positive.4[1,"pos.m"])
          start.snp         <- ibd.positive.4[1,"marker"]
          end.snp           <- ibd.positive.4[number.snps,"marker"]
          ibd.status        <- paste(unique(ibd.positive.4[,"viterbi"]), collapse=",")
          ibd.table <- rbind(ibd.table, cbind(fid1, iid1, fid2, iid2, chr, start.snp, end.snp, start.position.bp, end.position.bp, start.position.M, end.position.M, number.snps, length.bp, length.M, ibd.status))

        }
      }
    }
  }
  return(ibd.table)
}

