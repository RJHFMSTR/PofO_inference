#' XIBD Locus Proportion Plot
#'
#' Plot the proportion of pairs IBD for each SNP across the genome.
#' Annotation genes can be added to the plot and specific genes on interest can be highlighted.
#'
#' @param locus.proportions a data frame containing the proportion of pairs IBD at each SNP.
#' See \code{value} description in \code{\link{getLocusProportion}} for more details.
#' @param interval a vector of length 3 containing a chromosome, a start position (bp) and an end position (bp)
#' of an interval to plot. The default is \code{interval=NULL} which will plot the proportions over all chromosomes
#' in \code{locus.proportions}.
#' @param annotation.genes a data frame with at least 5 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{"numeric"} or \code{"integer"})
#' \item Gene name (type \code{"character"})
#' \item Start location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item End location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item Gene strand (+ or -)
#' }
#' \code{annotation.genes} should contain the following headers \code{chr, name, start, end} and \code{strand}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels.
#' @param highlight.genes a data frame with at least 4 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{"numeric"} or \code{"integer"})
#' \item Gene name (type \code{"character"})
#' \item Start location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item End location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' }
#' \code{highlight.genes} should contain the following headers \code{chr, name, start} and \code{end}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels.
#' @param add.rug Logical. Whether to include SNP positions as a rug in the figure. The default is \code{add.rug=FALSE}
#' @param plot.title a character string of a title to be added to the plot. The default is \code{plot.title=NULL} which does not add a title to the plot.
#' @import ggplot2
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
#'
#' # plot the proportion of pairs IBD
#' plotIBDproportions(locus.proportions = my_locus_prop,
#'                    interval = NULL,
#'                    annotation.genes = NULL,
#'                    highlight.genes = NULL,
#'                    add.rug = TRUE,
#'                    plot.title = NULL)
plotIBDproportions <- function(locus.proportions, interval = NULL, annotation.genes = NULL, highlight.genes = NULL, add.rug = TRUE, plot.title = NULL){

  # check locus matrix input
  if (ncol(locus.proportions) < 5)
    stop ("locus.proportions has incorrect format")
  colnames(locus.proportions)[1:4] <- c("chr","snp_id","pos_M","pos_bp")

  # if interval specified
  if (!is.null(interval) & length(interval) == 3) {
    chr   <- as.character(interval[1])
    start <- as.numeric(interval[2])
    stop  <- as.numeric(interval[3])
    if(!(chr %in% locus.proportions[,"chr"]))
      stop(paste0("chromosome ",chr," not in 'locus.proportions'\n"))
    if(start > stop)
      stop(paste("interval start=",interval[2]," is greater than interval end=",interval[3],sep=""))
    locus.interval <- locus.proportions[locus.proportions[,"chr"] == chr & locus.proportions[,"pos_bp"] >= start &
                                          locus.proportions[,"pos_bp"] <= stop,]
    if(nrow(locus.interval) == 0)
      stop("no SNPs in interval")
  } else {
    locus.interval <- locus.proportions
  }

  # check annotation.genes
  if (!is.null(annotation.genes)) {
    stopifnot(is.data.frame(annotation.genes))
    stopifnot(c("name","strand","chr","start","end") %in% colnames(annotation.genes))

    # subset genes if interval specified
    if (!is.null(interval)) {
      annotation.genes.0 <- annotation.genes[annotation.genes[,"chr"] == chr,]
      annotation.genes.1 <- NULL
      if(nrow(annotation.genes.0) > 0) {
        # function to find genes that overlap interval
        fun.overlap <- function(region.1, region.2) {
          a <- max(region.1[1], region.2[1])
          b <- min(region.1[2], region.2[2])
          return(c(a,b))
        }
        for(g in 1:nrow(annotation.genes.0)){
          gene.overlap <- fun.overlap(annotation.genes.0[g,c("start","end")],c(start,stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            annotation.genes.1 <- rbind(annotation.genes.1, annotation.genes.0[g,])
        }
      }
      annotation.genes.1 <- data.frame(annotation.genes.1)
    } else {
      annotation.genes.1 <- annotation.genes
    }
    if(nrow(annotation.genes.1) == 0) {
      cat("no annotation.genes in interval")
    } else {
      annotation.genes.1[,"start"] <- as.numeric(annotation.genes.1[,"start"])
      annotation.genes.1[,"end"] <- as.numeric(annotation.genes.1[,"end"])
    }
  }

  # check highlight.genes
  if (!is.null(highlight.genes)) {
    stopifnot(is.data.frame(highlight.genes))
    stopifnot(c("name","chr","start","end") %in% colnames(highlight.genes))

    # subset genes if interval specified
    if (!is.null(interval)) {
      highlight.genes.0 <- highlight.genes[highlight.genes[,"chr"] == chr,]
      highlight.genes.1 <- NULL
      if(nrow(highlight.genes.0) > 0) {
        # function to find genes that overlap interval
        fun.overlap <- function(region.1, region.2) {
          a <- max(region.1[1], region.2[1])
          b <- min(region.1[2], region.2[2])
          return(c(a,b))
        }
        for(g in 1:nrow(highlight.genes.0)){
          gene.overlap <- fun.overlap(highlight.genes.0[g,c("start","end")],c(start,stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            highlight.genes.1 <- rbind(highlight.genes.1, highlight.genes.0[g,])
        }
      }
      highlight.genes.1 <- data.frame(highlight.genes.1)
    } else {
      highlight.genes.1 <- highlight.genes
    }
    if(nrow(highlight.genes.1) == 0) {
      cat("no highlight.genes in interval")
    } else {
      highlight.genes.1[,"start"] <- as.numeric(highlight.genes.1[,"start"])
      highlight.genes.1[,"end"] <- as.numeric(highlight.genes.1[,"end"])
    }
  }

  # check add.rug
  stopifnot(is.logical(add.rug))

  # check title
  if (!is.null(plot.title)) {
    stopifnot(is.vector(plot.title))
    plot.title <- as.character(plot.title)
  }

  # chromosome names
  chromosomes <- as.character(unique(locus.interval[,"chr"]))

  # genome length
  genome.length <- 0
  for(chrom in chromosomes){
    positions.chrom <- locus.interval[locus.interval[,"chr"] == chrom,"pos_bp"]
    genome.length <- genome.length + (max(positions.chrom) - min(positions.chrom))
  }

  # define continious plot positions if interval not specified
  if (is.null(interval)) {
    chrstart <- 0
    chradd   <- 0
    labelpos <- NULL
    newpos   <- NULL
    for (i in 1:length(chromosomes)) {
      maxpos   <- max(locus.interval[locus.interval[,"chr"] == chromosomes[i],"pos_bp"])
      minpos   <- min(locus.interval[locus.interval[,"chr"] == chromosomes[i],"pos_bp"])
      newpos   <- c(newpos, locus.interval[locus.interval[,"chr"] == chromosomes[i],"pos_bp"] + chrstart)
      labelpos[i] <- (maxpos - minpos + 2*chrstart)/2
      chradd[i]   <- chrstart
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) != 0){
          annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"start"] <- annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"start"] + chrstart
          annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"end"] <- annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      if (!is.null(highlight.genes)) {
        if (nrow(highlight.genes.1) != 0){
          highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"start"] <- highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"start"] + chrstart
          highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"end"] <- highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      chrstart <- chrstart + maxpos
    }
    chradd <- chradd[-1]
  }

  # create dataframes for plotting
  if (is.null(interval)) {
    locus.df <- data.frame(pos=newpos,locus.interval[,5:ncol(locus.interval)])
    locus.df.1 <- locus.df
    locus.df.melt <- data.table::melt(locus.df.1,id="pos")
    colnames(locus.df.melt) <- c("pos","pairs","value")
    locus.df.melt[,"pos"] <- as.numeric(locus.df.melt[,"pos"])
    locus.df.melt[,"value"] <- as.numeric(locus.df.melt[,"value"])
  } else {
    locus.df <- data.frame(pos=locus.interval[,"pos_bp"],locus.interval[,5:ncol(locus.interval)])
    locus.df.melt <- data.table::melt(locus.df,id="pos")
    colnames(locus.df.melt) <- c("pos","pairs","value")
  }
  number.groups <- ncol(locus.df) - 1

  # plot:

  # setting up ggplot
  ggp <- ggplot()
  if (is.null(interval))
    ggp <- ggp + geom_vline(xintercept = chradd, colour = "gray87", linetype = "longdash")
  ggp <- ggp + geom_line(data = locus.df.melt, aes_string("pos", "value", col = "pairs"))
  ggp <- ggp + scale_colour_manual(values = getColourPaletteMajor(number.groups))
  ggp <- ggp + theme_bw()
  ggp <- ggp + ylab("Proportion of Pairs IBD")
  ggp <- ggp + theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   legend.title = element_blank(),
                   strip.background = element_rect(fill = "white",color = "white"),
                   strip.text.y = element_text(angle = 360))
  if (add.rug)
    ggp <- ggp + geom_rug(data=locus.df.melt, aes_string(x = "pos"), size = 0.1, colour = "gray30")
  if (!is.null(plot.title))
    ggp <- ggp + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
  if (number.groups == 1)
    ggp <- ggp + theme(legend.position = "none")

  # interval:
  if (is.null(interval)) {
    ggp <- ggp + xlab("Chromosome")
    #ggp <- ggp + geom_vline(xintercept = chradd, colour = "gray87", linetype = "longdash")
    ggp <- ggp + scale_x_continuous(breaks = labelpos, labels = chromosomes)
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  } else {
    ggp <- ggp + xlab(paste("Chromosome", chromosomes))
  }

  # annotation.genes:
  if (!is.null(annotation.genes)) {
    if(nrow(annotation.genes.1) != 0) {
      #min.y <- -0.05*max(locus.interval[,5:ncol(locus.interval)])
      #max.y <- -0.01*max(locus.interval[,5:ncol(locus.interval)])
      min.y <- -0.05*max(locus.df.melt[,"value"])
      max.y <- -0.01*max(locus.df.melt[,"value"])
      ggp <- ggp + geom_rect(data=annotation.genes.1, aes_string(xmin = "start", xmax = "end"), ymin = min.y, ymax = max.y, alpha = 0.9, fill = ifelse(annotation.genes.1[,"strand"]=="+","gold","red"))
      ggp <- ggp + ylim(min.y, max(locus.interval[,5:ncol(locus.interval)]))
    }
  }

  # highlight.genes:
  if (!is.null(highlight.genes)) {
    if(nrow(highlight.genes.1) != 0) {
      ggp <- ggp + geom_rect(data=highlight.genes.1, aes_string(xmin = "start", xmax = "end"), ymin = -Inf, ymax = Inf, fill = "gray40", alpha = 0.1)
      ggp <- ggp + geom_text(data=highlight.genes.1, aes_string(x = "start", y = 0, label = "name"), colour = "gray20", angle = 90, hjust = -0.1, vjust = -0.2, size = 3, alpha = 0.6)
      ggp <- ggp + geom_vline(data=highlight.genes.1, aes_string(xintercept = "start"), colour = "gray40", linetype = "solid", alpha = 0.1)
    }
  }

  print(ggp)
}
