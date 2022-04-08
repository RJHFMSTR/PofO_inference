//////////////////////////////////////////////////////////////////////
// rplot.cpp
// (c) 2010-2019 Wei-Min Chen
//
// This file is distributed as part of the KING source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile KING.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Oct 1, 2019

#include "rplot.h"
#include "MathMatrix.h"
#include "StringArray.h"

int CheckRout(const char *scriptfile);
void plotIBD1vsIBD2(const char *prefix, FILE *fp);

int plotAncestry(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_ancestryplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING Ancestry plot, by Zhennan Zhu and Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "library(e1071)\n");
   fprintf(fp, "prefix=\"%s\"\n", prefix);
   fprintf(fp, "pc <- read.table(paste0(prefix, \"pc.txt\"), header = TRUE)\n");
   fprintf(fp, "phe <- read.table(paste0(prefix, \"_popref.txt\"), header = TRUE)\n");
   fprintf(fp, "print(paste(\"Prepare the PC file and the reference file, starts at \",date()))\n");
   fprintf(fp, "pop <- phe[, c(\"IID\", \"Population\")]\n");
   fprintf(fp, "train.data <- pc[pc$AFF == 1, c(2, 7:16)]\n");
   fprintf(fp, "train.phe <- merge(train.data, pop, by = \"IID\")\n");
   fprintf(fp, "test.data <- pc[pc$AFF == 2, c(1, 2, 7:16)]\n");
   fprintf(fp, "train.x <- train.phe[, !colnames(train.phe) %%in%% c(\"Population\", \"IID\")]\n");
   fprintf(fp, "train.y <- train.phe[, \"Population\"]\n");
   fprintf(fp, "if(require(\"doParallel\", quietly=TRUE)){\n");
   fprintf(fp, "numCores <- detectCores()\n");
   fprintf(fp, "registerDoParallel(cores = min(round(numCores/2),41))\n");
   fprintf(fp, "tuneresults <- function(cost){\n");
   fprintf(fp, "tuneresult <- foreach(cost = cost, .combine = c) %%dopar%% {\n");
   fprintf(fp, "  set.seed(123)\n");
   fprintf(fp, "  mod = tune(svm, train.x, as.factor(train.y), kernel = \"linear\", cost = cost, probability = TRUE)\n");
   fprintf(fp, "  mod$performances[,c(\"error\")]\n");
   fprintf(fp, "}\n");
   fprintf(fp, "best.cost <- cost[which.min(tuneresult)]\n");
   fprintf(fp, "return(best.cost)}\n");
   fprintf(fp, "}else{\n");
   fprintf(fp, "  numCores <- 2\n");
   fprintf(fp, "  tuneresults <- function(cost) {\n");
   fprintf(fp, "    set.seed(123)\n");
   fprintf(fp, "    tune.mod = tune(svm, train.x, as.factor(train.y), kernel = \"linear\",\n");
   fprintf(fp, "      ranges=(list(cost=cost)), probability = TRUE)\n");
   fprintf(fp, "    return(tune.mod$best.parameters[1,1])}}\n");
   fprintf(fp, "print(paste(\"Assign\", min(round(numCores/2),41), \"cores for the grid search.\"))\n");
   fprintf(fp, "print(paste(\"Grid search with a wide range, starts at\", date()))\n");
   fprintf(fp, "best.cost <- tuneresults(2^(seq(-10, 10, by = 0.5)))\n");
   fprintf(fp, "print(paste(\"Grid search with a wide range, ends at\", date()))\n");
   fprintf(fp, "print(paste(\"The best cost is\", round(best.cost, 6), \"after the wide grid search\"))\n");
   fprintf(fp, "print(paste(\"Grid search with a small range, starts at\", date()))\n");
   fprintf(fp, "print(paste(\"The best cost is\", round(best.cost, 6), \"after the wide grid search\"))\n");
   fprintf(fp, "more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)\n");
   fprintf(fp, "best.cost <- tuneresults(more.cost)\n");
   fprintf(fp, "print(paste(\"Grid search with a small range, ends at\", date()))\n");
   fprintf(fp, "print(paste(\"The best cost is\", round(best.cost, 6), \"after the small grid search\"))\n");
   fprintf(fp, "set.seed(123)\n");
   fprintf(fp, "mymod <- svm(train.x, as.factor(train.y), cost = best.cost, kernel = \"linear\", probability=TRUE)\n");
   fprintf(fp, "print(paste(\"Predict ancestry information, start at\", date()))\n");
   fprintf(fp, "pred.pop <- predict(mymod, test.data[, !colnames(test.data) %%in%%c(\"FID\",\"IID\")], probability=TRUE)\n");
   fprintf(fp, "test.data$PRED <- pred.pop\n");
   fprintf(fp, "class.prob <- attr(pred.pop, \"probabilities\")\n");
   fprintf(fp, "print(paste(\"Prepare the summary file, starts at\", date()))\n");
   fprintf(fp, "orders <- t(apply(class.prob, 1, function(x) order(x,decreasing=T)))\n");
   fprintf(fp, "orders.class <- t(apply(orders, 1, function(x) colnames(class.prob)[x]))\n");
   fprintf(fp, "orders.probs <- t(sapply(1:nrow(class.prob), function(x) class.prob[x, orders[x,]]))\n");
   fprintf(fp, "check.cumsum <- t(apply(orders.probs, 1, cumsum))\n");
   fprintf(fp, "temp <- apply(check.cumsum, 1, function(x) which(x > 0.65)[1])\n");
   fprintf(fp, "pred_class <- sapply(1:length(temp), function (x) paste(orders.class[x, 1:as.numeric(temp[x])], collapse = \";\"))\n");
   fprintf(fp, "pred_prob <- sapply(1:length(temp), function (x) paste(round(orders.probs[x, 1:as.numeric(temp[x])], 3), collapse = \";\"))\n");
   fprintf(fp, "pred.out <- cbind(test.data[, c(\"FID\", \"IID\", \"PC1\", \"PC2\")], pred_class, pred_prob, orders.class[, 1], orders.class[, 2], round(orders.probs[, 1], 3), round(orders.probs[, 2], 3))\n");
   fprintf(fp, "colnames(pred.out)[5:10] <- c(\"Ancestry\", \"Pr_Anc\", \"Anc_1st\", \"Anc_2nd\", \"Pr_1st\", \"Pr_2nd\")\n");
   fprintf(fp, "print(paste(\"summary file is ready \", date()))\n");
   fprintf(fp, "write.table(pred.out, paste0(prefix, \"_InferredAncestry.txt\"), sep = \"\t\", quote = FALSE, row.names = FALSE)\n");
   fprintf(fp, "print(paste(\"Results are saved to\", paste0(prefix, \"_InferredAncestry.txt\")))\n");
   fprintf(fp, "print(\"Generate plots\")\n");
   fprintf(fp, "pred.out$Ancestry <- as.character(pred.out$Ancestry)\n");
   fprintf(fp, "pred.out$Ancestry[pred.out$Pr_1st <= 0.65] <- \">1 Pop\"\n");
   fprintf(fp, "Palette <- c(\"#1F78B4\",\"#33A02C\",\"#E31A1C\",\"#FF7F00\",\"#6A3D9A\",\"#B15928\",\"#A6CEE3\",\"#B2DF8A\",\"#FB9A99\",\"#FDBF6F\",\"#CAB2D6\",\"#FFFF99\",\"#999999\")\n");
   fprintf(fp, "train.groups <- unique(train.phe$Population)\n");
   fprintf(fp, "pred.colors <- rep(Palette[13], nrow(pred.out))\n");
   fprintf(fp, "for (i in 1:length(train.groups)){\n");
   fprintf(fp, "  pred.colors[pred.out$Ancestry==train.groups[i]] <- Palette[i]}\n");
   fprintf(fp, "train.colors <- rep(0, nrow(train.phe))\n");
   fprintf(fp, "for (i in 1:length(train.groups)){\n");
   fprintf(fp, "  train.colors[train.phe$Population==train.groups[i]] <- Palette[i]}\n");
   fprintf(fp, "x.adjust <- (max(train.phe$PC1, pred.out$PC1) - min(train.phe$PC1, pred.out$PC1))/10\n");
   fprintf(fp, "x.low <- min(train.phe$PC1, pred.out$PC1) - x.adjust\n");
   fprintf(fp, "x.high <- max(train.phe$PC1, pred.out$PC1) + x.adjust\n");
   fprintf(fp, "y.adjust <- (max(train.phe$PC2, pred.out$PC2) - min(train.phe$PC2, pred.out$PC2))/10\n");
   fprintf(fp, "y.low <- min(train.phe$PC2, pred.out$PC2) - y.adjust\n");
   fprintf(fp, "y.high <- max(train.phe$PC2, pred.out$PC2) + y.adjust\n");
   fprintf(fp, "postscript(paste0(prefix, \"_ancestryplot.ps\"), paper=\"letter\", horizontal=T)\n");
   fprintf(fp, "ncols <- min(3, ceiling(length(unique(pred.out$Ancestry))/2))\n");
   fprintf(fp, "if(!require(ggplot2, quietly=TRUE)) {\n");
   fprintf(fp, "  plot(pred.out$PC1, pred.out$PC2, col = pred.colors, xlab = \"PC1\", ylab = \"PC2\", xlim=c(x.low, x.high),\n");
   fprintf(fp, "    ylim=c(y.low, y.high), main=paste(\"Inferred Populations as Ancestry in\", prefix), pch = 16)\n");
   fprintf(fp, "  legend(\"topright\", legend = sort(unique(pred.out$Ancestry)), col = unique(pred.colors)[order(unique(pred.out$Ancestry))], pch = 16, cex = 1)\n");
   fprintf(fp, "par(mfrow=c(2,ncols))\n");
   fprintf(fp, "for (i in sort(unique(pred.out$Ancestry))) {\n");
   fprintf(fp, "  subdata <- subset(pred.out, Ancestry == i)\n");
   fprintf(fp, "  plot(subdata$PC1, subdata$PC2, col = unique(pred.colors)[unique(pred.out$Ancestry) == i],\n");
   fprintf(fp, "    xlim = c(x.low, x.high), ylim = c(y.low, y.high), xlab = \"PC1\", ylab = \"PC2\", main = paste0(i, \" (N=\", nrow(subdata), \")\"))}\n");
   fprintf(fp, "par(mfrow = c(1, 1))\n");
   fprintf(fp, "  plot(train.phe$PC1, train.phe$PC2, col = train.colors, xlim = c(x.low, x.high), ylim = c(y.low, y.high), xlab = \"PC1\", ylab = \"PC2\", main=\"Populations in Reference\", pch = 16)\n");
   fprintf(fp, "  legend(\"topright\", legend = sort(unique(train.phe$Population)),\n");
   fprintf(fp, "         col = unique(train.colors)[order(unique(train.phe$Population))], pch = 16, cex = 1)\n");
   fprintf(fp, "} else {\n");
   fprintf(fp, "  p <- ggplot(pred.out, aes(x = PC1, y = PC2))\n");
   fprintf(fp, "  p <- p + geom_point(aes(colour = factor(Ancestry, levels = sort(unique(Ancestry))))) +\n");
   fprintf(fp, "    xlim(x.low, x.high) + ylim(y.low, y.high) + labs(color = \"\") + \n");
   fprintf(fp, "    scale_colour_manual(values = unique(pred.colors)[order(unique(pred.out$Ancestry))]) +\n");
   fprintf(fp, "    ggtitle(paste(\"Inferred Populations as Ancestry in\", prefix))\n");
   fprintf(fp, "  print(p)\n");
   fprintf(fp, "  labels <- sapply(sort(unique(pred.out$Ancestry)), function(x) paste0(x, \" (N=\", sum(pred.out$Ancestry == x), \")\"))\n");
   fprintf(fp, "  p <- ggplot(pred.out, aes(x = PC1, y = PC2, colour = factor(Ancestry, levels = unique(Ancestry)))) +\n");
   fprintf(fp, "    scale_color_manual(values = unique(pred.colors)) + theme(legend.position = \"none\")\n");
   fprintf(fp, "  p <- p + geom_point() + xlim(x.low, x.high) + ylim(y.low, y.high) +\n");
   fprintf(fp, "    facet_wrap(~factor(Ancestry, levels = sort(unique(Ancestry)), labels=labels), ncol = min(3, ncols))\n");
   fprintf(fp, "  print(p)\n");
   fprintf(fp, "  p <- ggplot(train.phe, aes(x = PC1, y = PC2))\n");
   fprintf(fp, "  p <- p + geom_point(aes(colour = factor(Population, levels = sort(unique(Population))))) + xlim(x.low, x.high) + ylim(y.low, y.high)\n");
   fprintf(fp, "  p <- p + labs(color = \"\") + scale_colour_manual(values = unique(train.colors)[order(unique(train.phe$Population))]) + ggtitle(\"Populations in Reference\")\n");
   fprintf(fp, "  print(p)}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "print(paste0(prefix, \"_ancestryplot.ps is generated at \", date()))\n");
   fclose(fp);
   char command[256];
   if(rpath)
      sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
      sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("R plot for ancestry is not available for missing R libraries e1071.\n");
               printf("  Please rerun R code %s (or rerun KING) after e1071 is installed.\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_ancestryplot.ps", (const char*)prefix);
               system(command);
               printf("Ancestry populations are inferred as in %s_InferredAncestry.txt\n", (const char*)prefix);
               printf("Ancestry plots are generated in %s_ancestryplot.pdf\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
   return(error);
}

void plotROH(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_rohplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING ROH plot, by Zhennan Zhu, Jennifer N Nguyen, and Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "library(ggplot2)\n");
   fprintf(fp, "prefix=\"%s\"\n", prefix);
   fprintf(fp, "seg_name <- paste0(prefix, \".roh\")\n");
   fprintf(fp, "segments_name <- paste0(prefix, \".rohseg.gz\")\n");
   fprintf(fp, "all_seg_name <- paste0(prefix, \"allsegs.txt\")\n");
   fprintf(fp, "if( !(file.exists(seg_name) & file.exists(segments_name) & file.exists(all_seg_name)) ) stop(\"Missing RoH files\")\n");
   fprintf(fp, "all_seg <- read.table(all_seg_name, header = TRUE)\n");
   fprintf(fp, "all_seg <- subset(all_seg, select = c(Chr, StartMB, StopMB))\n");
   fprintf(fp, "segments <- read.table(segments_name, header = TRUE)\n");
   fprintf(fp, "segments <- subset(segments, select = c(FID, ID, Chr, StartMB, StopMB))\n");
   fprintf(fp, "roh <- read.table(seg_name, header = TRUE)\n");
   fprintf(fp, "roh_info <- roh[roh$F_ROH > 2^-4.5, c(\"FID\",\"ID\",\"F_ROH\")]\n");
   fprintf(fp, "roh_info$FID <- as.character(roh_info$FID)\n");
   fprintf(fp, "roh_info$ID <- as.character(roh_info$ID)\n");
   fprintf(fp, "postscript(paste0(prefix, \"_rohplot.ps\"), paper=\"letter\", horizontal = T)\n");
   fprintf(fp, "for(i in 1:nrow(roh_info)){\n");
   fprintf(fp, "  k <- subset(segments, FID == roh_info[i,1] & ID == roh_info[i,2])\n");
   fprintf(fp, "  if(nrow(k) > 0) {\n");
   fprintf(fp, "    theme_set(theme_bw(base_size = 18))\n");
   fprintf(fp, "    f_roh <- roh_info[i,\"F_ROH\"]\n");
   fprintf(fp, "    fid <- as.character(k[1,1])\n");
   fprintf(fp, "    id <- as.character(k[1,2])\n");
   fprintf(fp, "    g <- ggplot() +\n");
   fprintf(fp, "      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), fill = 'white', color = \"black\", size = 0.85) +\n");
   fprintf(fp, "      geom_rect(data = k, aes(xmin = StartMB, xmax = StopMB, ymin = 0, ymax = 0.9), fill = \"red\") +\n");
   fprintf(fp, "      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), color = \"black\", alpha = 0, size = 0.85) +\n");
   fprintf(fp, "      facet_grid(Chr ~ .) + scale_x_continuous(expand  = c(0, 0), limits = c(0, NA)) +\n");
   fprintf(fp, "      labs(x = \"Position (Mb)\", y = \"\", title = bquote(paste('Run of Homozygosity for ', .(id), ' from FAM ', .(fid), ' in ', .(prefix), ' (F'['ROH']*' = ', .(f_roh), ')'))) +\n");
   fprintf(fp, "      theme(legend.position = \"none\",\n");
   fprintf(fp, "        panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),\n");
   fprintf(fp, "        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),\n");
   fprintf(fp, "        axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title=element_text(size = 18))\n");
   fprintf(fp, "    print(g)}}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
      sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
      sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--roh is done but R plot is not available for missing R library ggplot2.\n");
               printf("  Please rerun R code %s (or rerun KING) after ggplot2 is installed.\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_rohplot.ps", (const char*)prefix);
               system(command);
               printf("ROH plots for inbred individuals are generated in %s_rohplot.pdf\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotMIerror(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_MIerrorplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING MI error plot, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "ped <- read.table(file=\"%ssplitped.txt\", stringsAsFactors=FALSE)[,c(1,2,5,6,7,8,9)]\n", (const char*)prefix);
   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0 | ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==2] <- 1\n");
   fprintf(fp, "colnames(ped) <- c(\"FID\", \"ID\", \"FA\", \"MO\", \"Sex\", \"Affected\", \"Status\")\n");
   fprintf(fp, "data <- read.table(\"%s.kin\",header = TRUE, stringsAsFactors = FALSE)\n", (const char*)prefix);
   fprintf(fp, "mi.err <- data[data[, \"Error\"]==1, \"FID\" ]\n");
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\", \"yellow\", NA)\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\", \"3rd\")\n");
   fprintf(fp, "data.fam <- merge(data, ped, by.x = c(\"FID\", \"ID1\"), by.y = c(\"FID\", \"ID\"))\n");
   fprintf(fp, "data.all <- merge(data.fam, ped, by.x = c(\"FID\", \"ID2\"), by.y = c(\"FID\", \"ID\"))[, c(\"FID\", \"ID1\", \"ID2\", \"Sex.x\", \"Sex.y\", \"InfType\", \"Error\")]\n");
   fprintf(fp, "data.all[!data.all$InfType %%in%% Inf.type, \"InfType\"] <- 6\n");
   fprintf(fp, "for (i in 1:5) data.all[data.all$InfType == Inf.type[i], \"InfType\"] <- i\n");
   fprintf(fp, "data.all[data.all$Error==1,\"Error\"] <- 5\n");
   fprintf(fp, "data.all[data.all$Error==0.5|data.all$Error==0,\"Error\"] <- 1\n");
   fprintf(fp, "postscript(\"%s_MIerrorplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", (const char*)prefix);
   fprintf(fp, "par(mfrow=c(1, 2))\n");
   fprintf(fp, "for (famid in unique(mi.err)){\n");
   fprintf(fp, "  fam <- ped[ped[, \"FID\"]==famid,]\n");
   fprintf(fp, "  if (all(fam[, c(\"FA\", \"MO\")]==0)){\n");
   fprintf(fp, "    g.empty <- make_empty_graph(n = nrow(fam), directed = FALSE)\n");
   fprintf(fp, "    plot(g.empty, vertex.size=27, vertex.color=NA, vertex.label.cex=.5, vertex.label.dist=1.6,\n");
   fprintf(fp, "         vertex.label.degree= pi/2, vertex.label.color=\"black\", vertex.label= fam[,\"ID\"],\n");
   fprintf(fp, "         edge.color=NA, layout=layout_with_fr(g.empty, grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "         vertex.shape=c(\"none\", \"square\", \"circle\")[1+fam[,\"Sex\"]])}else{\n");
   fprintf(fp, "  pedplot <- pedigree(id = fam$ID, dadid = fam$FA, momid = fam$MO, sex = as.numeric(fam$Sex),\n");
   fprintf(fp, "    affected = as.numeric(fam$Affected), status = as.numeric(fam$Status), famid = fam$FID, missid = 0)\n");
   fprintf(fp, "  plot(pedplot[toString(famid)], cex = 0.5, symbolsize = 2.8)}\n");
   fprintf(fp, "  fam.sub <- data.all[data.all$FID==famid,][, 2:7]\n");
   fprintf(fp, "  id <- unique(mapply(c, fam.sub[,c(1,3)], fam.sub[,c(2, 4)]))\n");
   fprintf(fp, "  g <- graph_from_data_frame(d=fam.sub, vertices=id[, 1], directed=FALSE)\n");
   fprintf(fp, "  plot(g, edge.width=fam.sub$Error, vertex.size=27, vertex.color=NA, vertex.label.cex=0.5,\n");
   fprintf(fp, "       edge.color=Inf.color[as.numeric(fam.sub$InfType)], layout=layout_with_fr(g, grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "       vertex.shape=c(\"none\", \"square\", \"circle\")[1+as.numeric(id[, 2])], margin=c(0.3,0,0,0))\n");
   fprintf(fp, "  legend(\"bottomright\", Inf.type, lty=1, col=Inf.color, text.col=Inf.color, cex=0.7, bty=\"n\")\n");
   fprintf(fp, "  mtext(paste(\"Documented versus Inferred Family\", famid, \"in %s\"), side = 3, line = -2, outer = TRUE)}\n", (const char*)prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--related is done but R plot is not available for missing R libraries kinship2/igraph.\n");
               printf("  Please rerun R code %s (or rerun KING) after kinship2/igraph is installed.\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_MIerrorplot.ps", (const char*)prefix);
               system(command);
               printf("MI error plots are generated in %s_MIerrorplot.pdf\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotUniqueFamily(const char *prefix, int degree, const char *analysis, const char *rpath)
{
   if(analysis != "related" && analysis != "ibdseg") return;
   String inputfile = prefix;
   if(analysis == "related") inputfile.Add(".kin0");
   else inputfile.Add(".seg");
   String namefile=prefix;
   if(analysis == "related") namefile.Add("_uniqfam0plot");
   else namefile.Add("_uniqfamplot");
   String scriptfile=namefile;
   scriptfile.Add(".R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --%s, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile, (const char*)analysis);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%s\", header=TRUE, stringsAsFactors=FALSE)[,c(\"ID1\", \"ID2\", \"InfType\")]\n", (const char*)inputfile);
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\")\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\")\n");
   fprintf(fp, "relatives <- data[data$InfType %%in%% Inf.type, ]\n");
   fprintf(fp, "for (i in 1:length(Inf.type)) relatives[relatives$InfType == Inf.type[i], \"InfType\"] <- i\n");
   fprintf(fp, "colnames(relatives)[3] <- \"color\"\n");
   fprintf(fp, "BuildPed <- function(v){\n");
   fprintf(fp, " buildType <- \"\"\n");
   fprintf(fp, " if(v[1]==3 && v[2]==3 && v[4]==2 && v[6]==1) buildType <- \"2_GG/HalfSibs\"\n");
   fprintf(fp, " else if(v[1]==3 && v[2]==3 && v[3]==1 && v[4]==2) buildType <- \"2_MZtwins+1_Parent/Child\"\n");
   fprintf(fp, " else if(v[1]==4 && v[2]==5 && v[3]==1 && v[4]==4) buildType <- \"2_Parents+2_MZtwins\"\n");
   fprintf(fp, " if(v[5]>0){\n");
   fprintf(fp, "  s <- floor(sqrt(v[5]*2))+1\n");
   fprintf(fp, "  if(s*(s-1)/2==v[5]){\n");
   fprintf(fp, "   if(v[2]==v[5] && v[4]==0) buildType <- paste0(s, \"_FullSiblings\")\n");
   fprintf(fp, "   else if(v[2]==v[5]+s && v[4]==s) buildType <- paste0(\"1_Parent+\", s, \"_FullSiblings\")\n");
   fprintf(fp, "   else if(v[2]==v[5]+s*2 && v[4]==s*2) buildType <- paste0(\"2_Parents+\", s, \"_FullSiblings\")\n");
   fprintf(fp, "   else buildType <- paste0(s, \"_FullSiblings+\", v[1]-s, \"_Relatives\")\n");
   fprintf(fp, "  }else if(v[3]==0 && v[4]==0 && v[6]==0) buildType <- \"Undetermined_FS_Only\"\n");
   fprintf(fp, "  else buildType <- \"Undetermined_FS\"\n");
   fprintf(fp, " }else if(v[4]==0 && v[6]==0){\n");
   fprintf(fp, "  buildType <- ifelse(v[3]==1, \"2_MZtwins\", \"Undetermined_MZ_Only\")\n");
   fprintf(fp, " }else if(v[3]==0 && v[4]==0){\n");
   fprintf(fp, "  buildType <- ifelse(v[6]==1, \"2_Second-Degree_Relatives\", \"Undetermined_2nd_Only\")\n");
   fprintf(fp, " }else if(v[3]==0 && v[6]==0){\n");
   fprintf(fp, "  if(v[4]==1) buildType <- \"1_Parent+1_Child\"\n");
   fprintf(fp, "  else if(v[4]==2) buildType <- \"2_Parents+1_Child(Trio)\"\n");
   fprintf(fp, "  else buildType <- \"Undetermined_PO_Only\"}\n");
   fprintf(fp, "if(buildType==\"\"){\n");
   fprintf(fp, " if(v[3]>0) buildType <- \"Undetermined_MZ\"\n");
   fprintf(fp, " else if(v[4]>0) buildType <- \"Undetermined_PO\"}\n");
   fprintf(fp, " return(buildType)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "g <- graph_from_data_frame(d = relatives, vertices = unique(c(relatives$ID1, relatives$ID2)), directed = FALSE)\n");
   fprintf(fp, "imc <- cluster_infomap(g)\n");
   fprintf(fp, "imc.membership <- membership(imc)\n");
   fprintf(fp, "community.info <- cbind(imc$membership,imc$names)\n");
   fprintf(fp, "colnames(community.info) <- c(\"membership\", \"ID\")\n");
   fprintf(fp, "community.rel <- merge(relatives, community.info, by.x = c(\"ID1\"), by.y = c(\"ID\"))\n");
   fprintf(fp, "community.rel$membership <- as.numeric(as.character(community.rel$membership))\n");
   fprintf(fp, "color_counts <- function(i) {\n");
   fprintf(fp, "  sub.colors <- community.rel[community.rel[, \"membership\"]==i, \"color\"]\n");
   fprintf(fp, "  return(sapply(1:4, function(y) sum(sub.colors==y)))}\n");
   fprintf(fp, "colors <- t(sapply(1:length(imc), color_counts))\n");
   fprintf(fp, "fam.info <- cbind(community.index=1:length(imc), sizes.node = sizes(imc), edge.colors=apply(colors, 1, sum), colors)\n");
   fprintf(fp, "fam.info <- as.data.frame(fam.info)\n");
   fprintf(fp, "all.counts <- aggregate(cbind(fam.info[0],counts=1),fam.info[,-1], length)\n");
   fprintf(fp, "fam.info.counts <- merge(fam.info, all.counts, by= c(\"sizes.node\",\"edge.colors\",\"V4\",\"V5\",\"V6\",\"V7\"))\n");
   fprintf(fp, "uniq.fam.info <- fam.info.counts[!duplicated(fam.info.counts[, c(1:6,8)]), ][,c(7, 1:6, 8)]\n");
   fprintf(fp, "uniq.fam.info <- uniq.fam.info[order(uniq.fam.info$counts,decreasing=TRUE),]\n");
   fprintf(fp, "all.counts <- uniq.fam.info[,\"counts\"]\n");
   fprintf(fp, "uniq.cluster <- induced_subgraph(g, (1:length(V(g)))[imc.membership %%in%% uniq.fam.info[, \"community.index\"]])\n");
   fprintf(fp, "all.names <- sapply(V(uniq.cluster)$name, function(x) membership(imc)[names(imc.membership)%%in%%x])\n");
   fprintf(fp, "all.builds <- apply((uniq.fam.info[, -1]),1,BuildPed)\n");
   fprintf(fp, "lo <- layout_(uniq.cluster, with_fr(), normalize())\n");
   fprintf(fp, "LocationForACluster<-function(x){\n");
   fprintf(fp, "  lo.local <- lo[all.names==x,]\n");
   fprintf(fp, "  return(c(min(as.numeric(lo.local[,1])), max(as.numeric(lo.local[,1])), max(as.numeric(lo.local[,2]))))}\n");
   fprintf(fp, "locations <- sapply(uniq.fam.info[, 1], LocationForACluster)\n");
   fprintf(fp, "postscript(\"%s.ps\", paper=\"letter\", horizontal=T)\n", (const char*)namefile);
   fprintf(fp, "plot(uniq.cluster, vertex.color=NA, vertex.size=1, vertex.label=NA, layout=lo, asp=0,\n");
   fprintf(fp, "  edge.color=Inf.color[as.numeric(E(uniq.cluster)$color)], main=\"All Unique Family Configurations in %s\")\n", prefix);
   fprintf(fp, "text((locations[1,]+locations[2,])/2, locations[3,]+0.04, all.counts)\n");
   fprintf(fp, "legend(\"bottomright\", Inf.type, lty = 1, col = Inf.color, text.col = Inf.color, cex = 0.7, bty = \"n\")\n");
   if(degree>1){
      fprintf(fp, "par(mfrow=c(2,3))\n");
      fprintf(fp, "is.built <- all.builds!=\"\"\n");
      fprintf(fp, "for(i in uniq.fam.info[, 1][is.built]){\n");
      fprintf(fp, "  index <- (1:length(uniq.fam.info[, 1]))[uniq.fam.info[, 1]==i]\n");
      fprintf(fp, "  g1 <- induced_subgraph(g, (1:length(V(g)))[imc.membership == i])\n");
      fprintf(fp, "  if(substr(all.builds[index],1,12)!=\"Undetermined\"&&all.counts[index]>1||all.counts[index]>10)\n");
      fprintf(fp, "  plot(g1, vertex.color=NA, vertex.size=3, vertex.label=NA, layout=layout_with_fr, asp=0,\n");
      fprintf(fp, "  edge.color=Inf.color[as.numeric(edge.attributes(g1)$color)], main=paste(all.counts[index], all.builds[index], \"Families\"))}\n");
   }
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   if(!CheckRout(scriptfile)){
      sprintf(command, "ps2pdf %s.ps", (const char*)namefile);
      system(command);
      printf("Unique family plot is generated in %s.pdf\n\n", (const char*)namefile);
   }
}

void plotDuplicate(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_duplicateplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --duplicate, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%s.con\", header=TRUE, stringsAsFactors=FALSE)[,c(2,4)]\n", prefix);
   fprintf(fp, "if(dim(data)[1]==0) q()\n");
   fprintf(fp, "postscript(\"%s_duplicateplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", prefix);
   fprintf(fp, "if(substr(data[1,1],1,4)==\"QRY_\"||substr(data[1,1],1,4)==\"REF_\"||substr(data[1,2],1,4)==\"QRY_\"||substr(data[1,2],1,4)==\"REF_\"){\n");
   fprintf(fp, " ordered <- rbind(data[substr(data[,1],1,3)==\"QRY\" & substr(data[,2],1,3)==\"REF\", c(1,2)],\n");
   fprintf(fp, "  data[substr(data[,1],1,3)==\"REF\" & substr(data[,2],1,3)==\"QRY\", c(2,1)])\n");
   fprintf(fp, " ordered <- cbind(substr(ordered[,1],5,1000),substr(ordered[,2],5,1000))\n");
   fprintf(fp, " mismatched <- ordered[ordered[,1]!=ordered[,2],]\n");
   fprintf(fp, " matched <- ordered[ordered[,1]==ordered[,2],]\n");
   fprintf(fp, " mismatched <- mismatched[!mismatched[,1] %%in%% matched[,1],]\n");
   fprintf(fp, " mismatched <- rbind(mismatched, matched[matched[,2] %%in%% mismatched[,2],])[,c(2,1)]\n");
   fprintf(fp, " if(dim(mismatched)[1]==0) q()\n");
   fprintf(fp, " g <- graph_from_data_frame(d=mismatched, vertices=unique(c(mismatched[,1],mismatched[,2])), directed=TRUE)\n");
   fprintf(fp, " plot(g, vertex.label.dist=0.3,vertex.label.degree=pi/2,vertex.size=1,vertex.label.cex=0.5, edge.arrow.size=0.5, edge.arrow.width=1,\n");
   fprintf(fp, "  layout=layout_with_fr, asp=0, main=paste(sum(mismatched[,1]!=mismatched[,2]), \"%s Samples Are Mismatched to Reference (QUERY<-REF)\"))\n", prefix);
   fprintf(fp, "}else{\n");
   fprintf(fp, "g <- graph_from_data_frame(d=data, vertices=unique(c(data$ID1,data$ID2)), directed=FALSE)\n");
   fprintf(fp, "if(dim(data)[1]<500) {plot(g, vertex.shape=\"none\", vertex.label.cex=0.5,\n");
   fprintf(fp, "  layout=layout_with_fr, asp=0, main=paste(dim(data)[1], \"Duplicate Pairs in %s\"))\n", prefix);
   fprintf(fp, "}else{g2 <- induced_subgraph(g, c((1:length(V(g)))[degree(g)>1],(1:length(V(g)))[degree(g)==1])[1:200])\n");
   fprintf(fp, "plot(g2, vertex.shape=\"none\", vertex.label.cex=0.5, layout=layout_with_fr, asp=0, main=\"200 Duplicates in %s\")}}\n", prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--duplicate is done but R plot is not available for missing R library igraph.\n");
               printf("  Please intall igraph and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--duplicate is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--duplicate is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_duplicateplot.ps", prefix);
               system(command);
               printf("Duplicate plot is generated in %s_duplicateplot.pdf\n\n", prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotCluster(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_clusterplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --cluster, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(igraph)\n");
   fprintf(fp, "data <- read.table(file=\"%scluster.kin\", header=TRUE, stringsAsFactors=FALSE)\n", (const char*)prefix);
   fprintf(fp, "postscript(\"%s_clusterplot.ps\", paper=\"letter\", horizontal=T, fonts=c(\"serif\", \"Palatino\"))\n", (const char*)prefix);
   fprintf(fp, "Inf.color <- c(\"purple\", \"red\", \"green\", \"blue\", \"yellow\")\n");
   fprintf(fp, "Inf.type <- c(\"Dup/MZ\", \"PO\", \"FS\", \"2nd\", \"3rd\")\n");
   fprintf(fp, "for(famid in unique(data$FID)){\n");
   fprintf(fp, "  fam <- data[data$FID==famid & data$InfType %%in%% Inf.type,][,c(2,3,4,5,15)]\n");
   fprintf(fp, "  id <- unique(mapply(c, fam[,c(1,3)], fam[,c(2,4)]))\n");
   fprintf(fp, "  for(i in 1:5) fam[fam$InfType==Inf.type[i],5] <- i\n");
   fprintf(fp, "  g <- graph_from_data_frame(d=fam, vertices=id[,1], directed=FALSE)\n");
   fprintf(fp, "  plot(g, edge.width=1.5, vertex.size=4, vertex.color=NA, vertex.label.cex=0.5,\n");
   fprintf(fp, "    edge.color=Inf.color[as.numeric(fam$InfType)], layout=layout_with_fr(g,grid=\"nogrid\"), asp=0,\n");
   fprintf(fp, "    vertex.shape=c(\"none\",\"square\",\"circle\")[1+as.numeric(id[,2])])\n");
   fprintf(fp, "  legend(\"bottomright\", Inf.type, lty=1, col=Inf.color, text.col=Inf.color, pt.cex=2, cex=0.8, bty=\"n\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--cluster is done but R plot is not available for missing R library igraph.\n");
               printf("  Please intall igraph and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--cluster is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--cluster is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_clusterplot.ps", (const char*)prefix);
               system(command);
               printf("Plots of newly clustered families are generated in %s_clusterplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotSplitped(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_pedplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING pedigree plot, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "ped <- read.table(file=\"%ssplitped.txt\", stringsAsFactors=FALSE)[,3:9]\n", (const char*)prefix);
   fprintf(fp, "postscript(\"%s_pedplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
//   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0] <- NA\n");
//   fprintf(fp, "ped$V8[ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==-9 | ped$V8==0 | ped$V8==1] <- 0\n");
   fprintf(fp, "ped$V8[ped$V8==2] <- 1\n");
   fprintf(fp, "pedAll <- pedigree(id = ped$V4, dadid = ped$V5, momid = ped$V6, sex = as.numeric(ped$V7), affected = as.numeric(ped$V8), status = as.numeric(ped$V9), famid = ped$V3, missid = 0)\n");
   fprintf(fp, "for(f in unique(ped$V3))\n");
   fprintf(fp, "  if(any(ped$V5[ped$V3 == f] != 0 | ped$V6[ped$V3 == f] != 0))\n");
   fprintf(fp, "    plot(pedAll[toString(f)], cex=0.5, symbolsize = 2.8)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("R plot is not available for missing R library kinship2.\n");
               printf("  Please rerun R code %s (or KING) after kinship2 is installed.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_pedplot.ps", (const char*)prefix);
               system(command);
               printf("Pedigree plots are generated in %s_pedplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotBuild(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_buildplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --build, by Wei-Min Chen and Zhennan Zhu\n", (const char*)scriptfile);
   fprintf(fp, "library(kinship2)\n");
   fprintf(fp, "ped <- read.table(file=\"%supdateparents.txt\", stringsAsFactors=FALSE)\n", (const char*)prefix);
   fprintf(fp, "ped$V6[ped$V6==-9 | ped$V6==0 | ped$V6==1] <- 0\n");
   fprintf(fp, "ped$V6[ped$V6==2] <- 1\n");
   fprintf(fp, "pedAll <- pedigree(id = ped$V2, dadid = ped$V3, momid = ped$V4, sex = as.numeric(ped$V5), affected = as.numeric(ped$V6), status = as.numeric(ped$V7),  famid = ped$V1, missid = 0)\n");
   fprintf(fp, "postscript(\"%s_buildplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "for(f in unique(ped$V1))\n");
   fprintf(fp, "  if(any(ped$V3[ped$V1 == f] != 0 | ped$V4[ped$V1 == f] != 0))\n");
   fprintf(fp, "    plot(pedAll[toString(f)], cex=0.5)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   switch(error){
      case 3:  printf("--build is done but R plot is not available for missing R library kinship2.\n");
               printf("  Please install kinship2 and run R code %s (or KING) again.\n\n", (const char*)scriptfile);
               break;
      case 2:  printf("--build is done but R code %s failed.\n\n", (const char*)scriptfile);
               break;
      case 1:  printf("--build is done but R code %s failed.\n", (const char*)scriptfile);
               printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
               break;
      case 0:  sprintf(command, "ps2pdf %s_buildplot.ps", (const char*)prefix);
               system(command);
               printf("Plots of newly reconstruction pedigrees are generated in %s_buildplot.pdf\n\n", (const char*)prefix);
               break;
      default: printf("Unexpected error. Please contact KING authors.\n");
   }
}

void plotHEreg(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_herplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --hereg, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_herplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.her\")) data <- read.table(file=\"%s.her\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "alltraits <- as.character(data$Trait)\n");
   fprintf(fp, "uniqtraits <- unique(alltraits)\n");
   fprintf(fp, "for(trait in uniqtraits){\n");
   fprintf(fp, "localdata <- data[alltraits==trait,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(0, min(max(LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = paste(\"Haseman-Elston Regression for Trait\" , trait))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "oldlocaldata <- localdata\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- oldlocaldata[oldlocaldata$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 2.2){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste0(i, \":\", localpos.max), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 2.2) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 2.2) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--hereg is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--hereg is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_herplot.ps", (const char*)prefix);
      system(command);
      printf("Haseman-Elston regression plots are generated in %s_herplot.pdf\n", (const char*)prefix);
   }
}

void plotPopROH(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_poprohplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --poproh, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_poprohplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.rohdiff\")) data <- read.table(file=\"%s.rohdiff\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "N <- max(data$Pop_Pos)\n");
   fprintf(fp, "label.pop <- 1:N\n");
   fprintf(fp, "for(p1 in 1:(N-1)){\n");
   fprintf(fp, "  for(p2 in (p1+1):N){\n");
   fprintf(fp, "localdata <- data[data$Pop_Neg == p1 & data$Pop_Pos == p2,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5, lwd = 1.5,\n");
   fprintf(fp, "  main = paste0(\"ROH difference in \", label.pop[p1], \" (-) and \", label.pop[p2], \" (+)\"),\n");
   fprintf(fp, "  ylim=c(max(min(LOD),-10), min(max(LOD),10)))\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,localdata$LOD[localdata$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(localdata$Pos[localdata$Chr==i][localdata$LOD[localdata$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD[localdata$Chr==i])\n");
   fprintf(fp, "if(minLOD <= -3.6){\n");
   fprintf(fp, "localpos <- median(localdata$Pos[localdata$Chr==i][localdata$LOD[localdata$Chr==i]==minLOD])\n");
   fprintf(fp, "text(base+localpos, max(minLOD,-10)+0.1, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--poproh is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--poproh is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_poprohplot.ps", (const char*)prefix);
      system(command);
      printf("ROH difference between populations is generated in %s_poprohplot.pdf\n", (const char*)prefix);
   }
}

void plotROHforQT(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_mthomoplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --mthomo, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_mthomoplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.mthomo\")) data <- read.table(file=\"%s.mthomo\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "alltraits <- as.character(data$Trait)\n");
   fprintf(fp, "uniqtraits <- unique(alltraits)\n");
   fprintf(fp, "for(trait in uniqtraits){\n");
   fprintf(fp, "localdata <- data[alltraits==trait,]\n");
   fprintf(fp, "Pos <- localdata$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[localdata$Chr==i] <- localdata$Pos[localdata$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos[localdata$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "LOD=localdata$LOD\n");
   fprintf(fp, "plot(Pos, LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(max(min(LOD),-10), min(max(LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = paste(\"Homozygosity Mapping of\" , trait))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "oldlocaldata <- localdata\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- oldlocaldata[oldlocaldata$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 2.2){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste0(i, \":\", localpos.max), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD)\n");
   fprintf(fp, "if(minLOD <= -2.2){\n");
   fprintf(fp, "localpos.min <- median(localdata$Pos[localdata$LOD==minLOD])\n");
   fprintf(fp, "text(base+localpos.min, max(minLOD,-10)+0.1, paste0(i, \":\", localpos.min), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 2.2) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(minLOD <= -2.2) localdata$LOD[localdata$Pos > localpos.min - 10 & localdata$Pos < localpos.min + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 2.2 && minLOD > -2.2) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = -2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--mthomo is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--mthomo is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_mthomoplot.ps", (const char*)prefix);
      system(command);
      printf("Homozygosity mapping plots are generated in %s_mthomoplot.pdf\n", (const char*)prefix);
   }
}

void plotROHmapping(const char *prefix, const char *stratName, int SEXCHR, const char *rpath)
{
   String postfix = stratName[0]=='\0'? "homomap": "homomapMH";
   String scriptfile=prefix;
   scriptfile.Add("_");
   scriptfile.Add(postfix);
   scriptfile.Add("plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --homomap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_%splot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix, (const char*)postfix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.%s\")) data <- read.table(file=\"%s.%s\", header=T)\n",
      (const char*)prefix, (const char*)postfix,
      (const char*)prefix, (const char*)postfix);
   fprintf(fp, "if(length(data$LOD)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LOD, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  ylim=c(max(min(data$LOD),-5), min(max(data$LOD),10)),\n");
   fprintf(fp, "  lwd = 1.5, main = \"Homozygosity Mapping in %s\")\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "localdata <- data[data$Chr==i,]\n");
   fprintf(fp, "repeat{\n");
   fprintf(fp, "maxLOD <- max(0,localdata$LOD)\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos.max <- median(localdata$Pos[localdata$LOD==maxLOD])\n");
   fprintf(fp, "text(base+localpos.max, min(maxLOD,10)-0.1, paste0(i, \":\", localpos.max), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "minLOD <- min(0,localdata$LOD)\n");
   fprintf(fp, "if(minLOD <= -3.6){\n");
   fprintf(fp, "localpos.min <- median(localdata$Pos[localdata$LOD==minLOD])\n");
   fprintf(fp, "text(base+localpos.min, max(minLOD,-5)+0.1, paste0(i, \":\", localpos.min), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(maxLOD >= 3.6) localdata$LOD[localdata$Pos > localpos.max - 10 & localdata$Pos < localpos.max + 10] <- 0;\n");
   fprintf(fp, "if(minLOD <= -3.6) localdata$LOD[localdata$Pos > localpos.min - 10 & localdata$Pos < localpos.min + 10] <- 0;\n");
   fprintf(fp, "if(maxLOD < 3.6 && minLOD > -3.6) break\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, localdata$Pos)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "abline(h = -3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--homomap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--homomap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_%splot.ps",
         (const char*)prefix, (const char*)postfix);
      system(command);
      printf("Homozygosity mapping plot is generated in %s_%splot.pdf\n",
         (const char*)prefix, (const char*)postfix);
   }
}

void plotPopDist(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_popdistplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --popdist, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_popdistplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.dst\")) data <- read.table(file=\"%s.dst\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$DistIBD)>0){\n");
   fprintf(fp, "N <- sqrt(2*length(data$DistIBD)+0.25)-0.5\n");
   fprintf(fp, "M <- matrix(0,N,N)\n");
   fprintf(fp, "M[upper.tri(M, diag=TRUE)] <- data$DistIBD\n");
   fprintf(fp, "M[lower.tri(M, diag=TRUE)] <- data$DistIBD\n");
   fprintf(fp, "plot(hclust(as.dist(M)),main=\"Clustering of %s Populations By Average IBD\",xlab=\"\",ylab=\"1 - IBD Proportion\")\n",
      (const char*)prefix);
   fprintf(fp, "mds <- cmdscale(M)\n");
   fprintf(fp, "plot(mds[,1],mds[,2],type=\"n\",xlab=\"\",ylab=\"\", axes=FALSE,main=\"MDS of %s Populations Using Average IBD\")\n",
      (const char*)prefix);
   fprintf(fp, "text(mds[,1],mds[,2],1:N)\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--popdist is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--popdist is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_popdistplot.ps", (const char*)prefix);
      system(command);
      printf("Population distance plot is generated in %s_popdistplot.pdf\n", (const char*)prefix);
   }
}

void plotIBDmapping(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibdmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ibdmap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibdmapplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.ibdmap\")) data <- read.table(file=\"%s.ibdmap\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$P)>0){\n");
   fprintf(fp, "Pos <- data$PosMb\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$PosMb[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0, data$PosMb[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "logP <- rep(10, length(data$P))\n");
   fprintf(fp, "logP[data$P>0] <- -log(data$P[data$P>0])/log(10)\n");
   fprintf(fp, "plot(Pos, logP, type=\"l\", xlab=\"Chromosome\", ylab = expression(paste(\"-\", log[10],\"P\",sep=\"\")), xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.2, cex.axis=1.2, cex.main=1.5, cex.lab = 1.2, cex.sub = 1.2, ylim=c(0,max(c(logP[logP<10]),4)),\n");
   fprintf(fp, "  lwd = 1.2, main = \"Permutation P Values of IBD Mapping in %s\")\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "minP <- min(1,data$P[data$Chr==i])\n");
   fprintf(fp, "if(minP < 0.0001){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$P[data$Chr==i]==minP])\n");
   fprintf(fp, "text(base+localpos, -log(minP, 10)-0.1, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 4, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ibdmap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdmap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ibdmapplot.ps", (const char*)prefix);
      system(command);
      printf("IBD mapping plot is generated in %s_ibdmapplot.pdf\n", (const char*)prefix);
   }
}

/*
void plotAncestry(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ancestryplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ancestry, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ancestryplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.anc\")) data <- read.table(file=\"%s.anc\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$Admix)>0){\n");
   fprintf(fp, "  plot(data$Admix, data$Ancestry, xlab=\"Admixed Proportion\", ylab = \"Ancestry\",\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"Ancestry in %s\")\n", (const char*)prefix);
   fprintf(fp, "  abline(a=0, b=0.5, col=\"red\")\n");
   fprintf(fp, "  abline(a=1, b=-0.5, col=\"red\")\n");
   fprintf(fp, "  y <- seq(0,1,0.001)\n");
   fprintf(fp, "  x <- y*(1-y)*2\n");
   fprintf(fp, "  lines(x, y, col=\"red\")\n");
   fprintf(fp, "  plot(data$Anc_P1, data$Anc_P2, xlab=\"Ancestry of Parent 1\", ylab = \"Ancestry of Parent 2\",\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"Parental Ancestry in %s\")\n", (const char*)prefix);
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ancestry is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ancestry is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ancestryplot.ps", (const char*)prefix);
      system(command);
      printf("Ancestry plot is generated in %s_ancestryplot.pdf\n", (const char*)prefix);
   }
}
*/

void plotNPL(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_nplplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --npl, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_nplplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.npl\")) data <- read.table(file=\"%s.npl\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$LODwDSP)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0,data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LODwDSP, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"ASP/DSP NPL Scan in %s\", ylim=c(0,min(max(data$LODwDSP),10)))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,data$LODwDSP[data$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$LODwDSP[data$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(0,data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos, data$LOD_ASP, type=\"l\", xlab=\"Chromosome\", ylab = \"LOD Score\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5,\n");
   fprintf(fp, "  lwd = 1.5, main = \"ASP NPL Scan in %s\", ylim=c(0,min(max(data$LOD_ASP),10)))\n", (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxLOD <- max(0,data$LOD_ASP[data$Chr==i])\n");
   fprintf(fp, "if(maxLOD >= 3.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$LOD_ASP[data$Chr==i]==maxLOD])\n");
   fprintf(fp, "text(base+localpos, min(maxLOD,10)-0.1, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 3.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "abline(h = 2.2, lty = 3, col=\"green\")\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--npl is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--npl is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_nplplot.ps", (const char*)prefix);
      system(command);
      printf("ASP NPL plot is generated in %s_nplplot.pdf\n", (const char*)prefix);
   }
}

void plotAUCmapping(const char *prefix, int SEXCHR, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_aucmapplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --aucmap, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_aucmapplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.aucmap\")) data <- read.table(file=\"%s.aucmap\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "if(length(data$AUC)>0){\n");
   fprintf(fp, "Pos <- data$Pos\n");
   fprintf(fp, "baseStop <- c()\n");
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "  Pos[data$Chr==i] <- data$Pos[data$Chr==i] + base\n");
   fprintf(fp, "  base <- base + max(data$Pos[data$Chr==i])\n");
   fprintf(fp, "  baseStop <- c(baseStop, base)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "baseStart <- c(0, baseStop[-%d])\n", SEXCHR-1);
   fprintf(fp, "baseMiddle <- (baseStart + baseStop)/2\n");
   fprintf(fp, "plot(Pos[data$Success>0.9], data$AUC[data$Success>0.9], type=\"l\", xlab=\"Chromosome\", ylab = \"AUC\", xaxt = 'n',\n");
   fprintf(fp, "  cex = 1.5, cex.axis=1.5, cex.main=1.5, cex.lab = 1.5, cex.sub = 1.5, ylim=c(0.5,max(data$AUC[data$Success>0.9])),\n");
   fprintf(fp, "  lwd = 1.5, main = \"Risk Prediction Using Position-Specific IBD Relatives in %s\")\n",
      (const char*)prefix);
   fprintf(fp, "base <- 0\n");
   fprintf(fp, "for(i in 1:%d){\n", SEXCHR-1);
   fprintf(fp, "maxAUC <- max(0,data$AUC[data$Chr==i])\n");
   fprintf(fp, "if(maxAUC > 0.6){\n");
   fprintf(fp, "localpos <- median(data$Pos[data$Chr==i][data$AUC[data$Chr==i]==maxAUC])\n");
   fprintf(fp, "text(base+localpos, maxAUC-0.01, paste0(i, \":\", localpos), col=\"red\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "  base <- base + max(0, data$Pos[data$Chr==i])\n");
   fprintf(fp, "}\n");
   fprintf(fp, "axis(1, labels=FALSE, line = 0, at = c(0,baseStop))\n");
   fprintf(fp, "axis(1, labels=c(1:%d), line = -0.5, lty=\"blank\", at = baseMiddle)\n", SEXCHR-1);
   fprintf(fp, "abline(h = 0.6, lty = 2, col=\"red\")\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "}\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--aucmap is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--aucmap is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_aucmapplot.ps", (const char*)prefix);
      system(command);
      printf("AUC mapping plot is generated in %s_aucmapplot.pdf\n", (const char*)prefix);
   }
}

void plotRelationship(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_relplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --related, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_relplot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.kin\")) data <- read.table(file=\"%s.kin\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   fprintf(fp, "allcolors <- c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\")\n");
   fprintf(fp, "allrelatives <- c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"3rd-Degree\", \"More Distant\", \"Unrelated\")\n");
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   // Page 1 plot: IBD1 vs IBD2
   fprintf(fp, "allpair <- data$PropIBD>0 | data$Kinship>0.04419\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi>0.08839 & data$Phi<=0.17678\n");
   fprintf(fp, "d3 <- data$Phi>0.04419 & data$Phi<=0.08839\n");
   fprintf(fp, "dO <- data$Phi>0 & data$Phi<=0.04419\n");
   fprintf(fp, "dU <- data$Phi==0 & allpair\n");
   fprintf(fp, "plot(data$IBD1Seg[dU], data$IBD2Seg[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$IBD1Seg[allpair]), max(data$IBD1Seg[allpair])),\n");
   fprintf(fp, "ylim=c(min(data$IBD2Seg[allpair]), max(data$IBD2Seg[allpair])),\n");
   fprintf(fp, "main = \"IBD Segments In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Length Proportion of IBD1 Segments (\", pi[1], \")\", sep=\"\")),\n");
   fprintf(fp, "ylab=expression(paste(\"Length Proportion of IBD2 Segments (\", pi[2], \")\", sep=\"\")))\n");
   fprintf(fp, "points(data$IBD1Seg[d0], data$IBD2Seg[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO], data$IBD2Seg[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS], data$IBD2Seg[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d2], data$IBD2Seg[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$IBD1Seg[d3], data$IBD2Seg[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$IBD1Seg[dO], data$IBD2Seg[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$IBD1Seg[dU & data$PropIBD>0.08838835], data$IBD2Seg[dU & data$PropIBD>0.08838835], col=\"black\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS & data$IBD2Seg<0.08], data$IBD2Seg[d1.FS & data$IBD2Seg<0.08], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO & data$IBD1Seg+data$IBD2Seg<0.9], data$IBD2Seg[d1.PO & data$IBD1Seg+data$IBD2Seg<0.9], col=\"red\")\n");
   fprintf(fp, "abline(h = 0.08, col = \"green\", lty = 3, lwd = 2)\n");   // FS && MZ
   fprintf(fp, "abline(a = 0.96, b = -1, col = \"red\", lty = 3, lwd = 2)\n");   // PO
   fprintf(fp, "abline(a = 0.3535534, b = -0.5, col = \"green\", lty = 3, lwd = 2)\n");// FS
   fprintf(fp, "abline(a = 0.1767767, b = -0.5, col = \"blue\", lty = 3, lwd = 2)\n");// 2nd/3rd
   fprintf(fp, "abline(a = 0.08838835, b = -0.5, col = \"magenta\", lty = 3, lwd = 2)\n");// 3rd/4th
   fprintf(fp, "abline(a = 0.04419, b = -0.5, col = \"gold\", lty = 3, lwd = 2)\n");// 4th/UN
   fprintf(fp, "legend(\"topright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");

   // Page 2 plot: Kinship vs. PropIBD
   fprintf(fp, "allpair <- data$PropIBD>0 | data$Kinship>0.04419\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi>0.08839 & data$Phi<=0.17678\n");
   fprintf(fp, "d3 <- data$Phi>0.04419 & data$Phi<=0.08839\n");
   fprintf(fp, "dO <- data$Phi>0 & data$Phi<=0.04419\n");
   fprintf(fp, "dU <- data$Phi==0 & allpair\n");
   fprintf(fp, "plot(data$PropIBD[dU], data$Kinship[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$PropIBD[allpair]), max(data$PropIBD[allpair])),\n");
   fprintf(fp, "ylim=c(min(data$Kinship[allpair]), max(data$Kinship[allpair])),\n");
   fprintf(fp, "main = paste(\"Kinship vs Proportion IBD (Corr=\", round(cor(data$Kinship[data$PropIBD>0], data$PropIBD[data$PropIBD>0]),digit=3),\") in %s Families\",sep=\"\"),\n",
      (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Proportion of Genomes IBD (\", pi,\"=\",pi[2]+pi[1]/2,\")\",sep=\"\")),\n");
   fprintf(fp, "ylab = expression(paste(\"Estimated Kinship Coefficient (\", phi, \")\",sep=\"\")))\n");
   fprintf(fp, "points(data$PropIBD[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$PropIBD[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$PropIBD[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$PropIBD[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$PropIBD[d3], data$Kinship[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$PropIBD[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$PropIBD[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.35355, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.17678, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.08838, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.04419, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.02210, col = \"gold\", lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.5, lty = 1)\n");
   fprintf(fp, "abline(a = 0, b = 0.7071068, lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.3535534, lty = 3)\n");
   fprintf(fp, "abline(v = 0.70711, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.35355, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.17678, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.08838, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.04419, col = \"gold\", lty = 3)\n");
   fprintf(fp, "text(x=0.35355, y=0.35355, \"1st\", adj=c(0,1), col=\"green\")\n");
   fprintf(fp, "text(x=0.17678, y=0.17678, \"2nd\", adj=c(0,1), col=\"blue\")\n");
   fprintf(fp, "text(x=0.08839, y=0.08839, \"3rd\", adj=c(0,1), col=\"magenta\")\n");
   fprintf(fp, "text(x=0.04419, y=0.04419, \"4th\", adj=c(0,1), col=\"gold\")\n");
   fprintf(fp, "legend(\"bottomright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");

   fprintf(fp, "}else if(length(data$Kinship)>0){\n");
   // In absence of IBDSeg, plot 1: Kinship vs HomIBS0
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi==0.125\n");
   fprintf(fp, "dU <- data$Phi==0\n");
   fprintf(fp, "dO <- !d0 & !d1.PO & !d1.FS & !d2 & !dU\n");
   fprintf(fp, "plot(data$HomIBS0[dU], data$Kinship[dU], type=\"p\", col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim=c(min(data$HomIBS0), max(data$HomIBS0)),\n");
   fprintf(fp, "ylim=c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = \"Kinship vs IBS0 in %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of Zero IBS In Minor Homozygote Pairs\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$HomIBS0[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$HomIBS0[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$HomIBS0[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$HomIBS0[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$HomIBS0[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$HomIBS0[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"topright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");
   // in absence of IBDSeg, plot 2: Kinship vs. HetConc
   fprintf(fp, "plot(data$HetConc[dU], data$Kinship[dU], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim = c(min(data$HetConc), max(data$HetConc[data$HetConc<0.8])),\n");
   fprintf(fp, "ylim = c(min(data$Kinship), max(data$Kinship[data$HetConc<0.8])),\n");
   fprintf(fp, "main = \"Kinship vs Heterozygote Concordance In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Heterozygote Concordance Rate\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$HetConc[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$HetConc[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$HetConc[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$HetConc[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$HetConc[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"bottomright\", allrelatives, col=allcolors, text.col=allcolors, pch=19, cex=1.2)\n");
   fprintf(fp, "if(sum(d1.FS)>20){\n");
   fprintf(fp, "y.cut <- 2^-2.5\n");
   fprintf(fp, "x.FS <- quantile(data$HetConc[d1.FS], probs=c(0.1, 0.9))\n");
   fprintf(fp, "y.FS <- quantile(data$Kinship[d1.FS], probs=c(0.1, 0.9))\n");
   fprintf(fp, "slope.FS <- (y.FS[2]-y.FS[1]) / (x.FS[2]-x.FS[1])\n");
   fprintf(fp, "a.FS <- y.FS[1]-x.FS[1]*slope.FS\n");
   fprintf(fp, "abline(a=a.FS, b=slope.FS, col=\"green\")\n");
   fprintf(fp, "x.cut.FS <- (y.cut - a.FS)/slope.FS\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(sum(d2)>20){\n");
   fprintf(fp, "x.d2 <- quantile(data$HetConc[d2], probs=c(0.9, 0.1))\n");
   fprintf(fp, "y.d2 <- quantile(data$Kinship[d2], probs=c(0.9, 0.1))\n");
   fprintf(fp, "slope.d2 <- (y.d2[2]-y.d2[1]) / (x.d2[2]-x.d2[1])\n");
   fprintf(fp, "a.d2 <- y.d2[1]-x.d2[1]*slope.d2\n");
   fprintf(fp, "abline(a=a.d2, b=slope.d2, col=\"blue\")\n");
   fprintf(fp, "x.cut.d2 <- (y.cut - a.d2)/slope.d2\n");
   fprintf(fp, "x.cut <- sqrt(x.cut.d2 * x.cut.FS)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "if(sum(d1.FS)>20 & sum(d2)>20){\n");
   fprintf(fp, "print(c(x.cut.d2, x.cut.FS, x.cut))\n");
   fprintf(fp, "abline(v=x.cut, col=\"purple\", lty=3)\n");
   fprintf(fp, "text(x=x.cut, y=0.3, labels=round(x.cut,digit=3),col=\"purple\")\n");
   fprintf(fp, "points(data$HetConc[d1.FS & data$HetConc < x.cut], data$Kinship[d1.FS & data$HetConc < x.cut], col=\"green\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "}\n");
   fprintf(fp, "\nif(!file.exists(\"%s.kin0\")) quit()\n", (const char*)prefix);
   fprintf(fp, "data <- read.table(file = \"%s.kin0\", header=T)\n", (const char*)prefix);
   plotIBD1vsIBD2(prefix, fp);

   // Page 2 plot: Kinship vs. PropIBD
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   fprintf(fp, "d0 <- data$IBD2Seg>0.7\n");
   fprintf(fp, "d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)\n");
   fprintf(fp, "d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08\n");
   fprintf(fp, "d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)\n");
   fprintf(fp, "d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678\n");
   fprintf(fp, "d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839\n");
   fprintf(fp, "dU <- data$PropIBD<=0.04419\n");
   fprintf(fp, "plot(data$PropIBD[d1.FS], data$Kinship[d1.FS], type=\"p\", col=\"green\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$PropIBD), max(data$PropIBD)),\n");
   fprintf(fp, "ylim=c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = paste(\"Kinship vs Proportion IBD (Corr=\", round(cor(data$Kinship[data$PropIBD>0], data$PropIBD[data$PropIBD>0]),digit=3),\") in Inferred %s Relatives\",sep=\"\"),\n",
      (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Proportion of Genomes IBD (\", pi,\"=\",pi[2]+pi[1]/2,\")\",sep=\"\")), \n");
   fprintf(fp, "ylab = expression(paste(\"Estimated Kinship Coefficient (\", phi, \")\",sep=\"\")))\n");
   fprintf(fp, "points(data$PropIBD[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$PropIBD[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$PropIBD[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$PropIBD[d3], data$Kinship[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$PropIBD[d4], data$Kinship[d4], col=\"gold\")\n");
   fprintf(fp, "points(data$PropIBD[dU], data$Kinship[dU], col=\"black\")\n");
   fprintf(fp, "abline(h = 0.35355, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.17678, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.08839, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.04419, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.02210, col = \"gold\", lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.5, lty = 1)\n");
   fprintf(fp, "abline(a = 0, b = 0.7071068, lty = 3)\n");
   fprintf(fp, "abline(a = 0, b = 0.3535534, lty = 3)\n");
   fprintf(fp, "abline(v = 0.70711, col = \"purple\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.35355, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.17678, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.08839, col = \"magenta\", lty = 3)\n");
   fprintf(fp, "abline(v = 0.04419, col = \"gold\", lty = 3)\n");
   fprintf(fp, "text(x=0.35355, y=0.35355, \"1st\", adj=c(0,1), col=\"green\")\n");
   fprintf(fp, "text(x=0.17678, y=0.17678, \"2nd\", adj=c(0,1), col=\"blue\")\n");
   fprintf(fp, "text(x=0.08839, y=0.08839, \"3rd\", adj=c(0,1), col=\"magenta\")\n");
   fprintf(fp, "text(x=0.04419, y=0.04419, \"4th\", adj=c(0,1), col=\"gold\")\n");
   fprintf(fp, "legend(\"bottomright\", c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "col=allcolors, text.col = allcolors, pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--related is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--related is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_relplot.ps", (const char*)prefix);
      system(command);
      printf("Relationship plot is generated in %s_relplot.pdf\n", (const char*)prefix);
   }
}

void plotIBDSeg(const char *prefix, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd1vsibd2.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --ibdseg, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd1vsibd2.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data<-c()\n");
   fprintf(fp, "if(file.exists(\"%s.seg\")) data <- read.table(file=\"%s.seg\", header=T)\n",
      (const char*)prefix,(const char*)prefix);
   plotIBD1vsIBD2(prefix, fp);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--ibdseg is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdseg is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_ibd1vsibd2.ps", (const char*)prefix);
      system(command);
      printf("IBD1 vs IBD2 plot is generated in %s_ibd1vsibd2.pdf\n", (const char*)prefix);
   }
}

void plotGenderError(const char *prefix, IntArray & plotx, Vector & ploty, IntArray & plotz, double xHeterozygosity, int gendererrorCount, const char *rpath)
{
   String genderplotdata=prefix;
   genderplotdata.Add("_gender_autodata.txt");
   FILE *fp = fopen(genderplotdata, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)genderplotdata);
   fprintf(fp, "YCount\txHeterozygosity\tSEX\n");
   for(int i = 0; i < plotx.Length(); i++)
      fprintf(fp, "%d\t%.5lf\t%d\n", plotx[i], ploty[i], plotz[i]);
   fclose(fp);
   String genderplot=prefix;
   genderplot.Add("_gender_autoplot.ps");
   String scriptfile=prefix;
   scriptfile.Add("_gender_autoplot.R");
   fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --autoQC, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "data <- read.table(file=\"%s\", header=T)\n", (const char*)genderplotdata);
   fprintf(fp, "postscript(\"%s\", paper=\"letter\", horizontal=T)\n", (const char*)genderplot);
   fprintf(fp, "isFemale<-data$SEX==2\n");
   fprintf(fp, "isMale<-data$SEX==1\n");
   fprintf(fp, "isUnknown<-data$SEX==0\n");
   fprintf(fp, "cutoff<-max(data$YCount)/2\n");
   fprintf(fp, "plot(data$YCount[isFemale],data$xHeterozygosity[isFemale], type=\"p\",\n");
   fprintf(fp, "  col=\"red\", xlim=c(0,max(data$YCount)), ylim=c(0,max(data$xHeterozygosity)),\n");
   if(gendererrorCount){
      fprintf(fp, "  main=\"Gender Checking in %s Samples (%d Samples Mislabeled)\", \n",
         (const char*)prefix, gendererrorCount);
      fprintf(fp, "  xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }else{
      fprintf(fp, "  main=\"Gender Checking in %s Samples\", \n", (const char*)prefix);
      fprintf(fp, "  xlab=\"# Y-Chr SNPs\", ylab=\"X-Chr Heterozygosity\")\n");
   }
   fprintf(fp, "points(data$YCount[isMale], data$xHeterozygosity[isMale], col=\"blue\")\n");
   fprintf(fp, "points(data$YCount[isUnknown], data$xHeterozygosity[isUnknown], col=\"black\")\n");
   fprintf(fp, "points(data$YCount[isFemale&data$YCount>cutoff], data$xHeterozygosity[isFemale&data$YCount>cutoff], col=\"red\")\n");
   fprintf(fp, "points(data$YCount[isMale&data$YCount<cutoff], data$xHeterozygosity[isMale&data$YCount<cutoff], col=\"blue\")\n");
   fprintf(fp, "abline(v=cutoff, col=\"purple\", lty=2)\n");
   fprintf(fp, "abline(v=cutoff*2/3, col=\"purple\")\n");
   fprintf(fp, "abline(v=cutoff*4/3, col=\"purple\")\n");
   fprintf(fp, "abline(h=%.4lf, col=\"purple\")\n", xHeterozygosity);
   fprintf(fp, "legend(\"topright\", c(\"Female\", \"Male\", \"Unknown\"), col=c(\"red\", \"blue\", \"black\"),\n");
   fprintf(fp, "  text.col = c(\"red\", \"blue\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--autoQC is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--autoQC is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s", (const char*)genderplot);
      system(command);
      printf("Gender plot are generated in %s_gender_autoplot.pdf\n", (const char*)prefix);
   }
}

void plotPopStructure(const char *prefix, int projectFlag, const char *rpath)
{
   String scriptfile=prefix;
   scriptfile.Add("_pcplot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for KING --mds or --pca, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_pcplot.ps\", paper=\"letter\", horizontal=T)\n",(const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%spc.txt\", header=T)\n", (const char*)prefix);
   fprintf(fp, "plot(data$PC1, data$PC2, type=\"p\", xlab=\"PC1\", ylab=\"PC2\", main = \"Population Structure in %s\")\n",(const char*)prefix);
   if(projectFlag){
      fprintf(fp, "points(data[data[,6]==2,7], data[data[,6]==2,8], col = \"red\")\n");
      fprintf(fp, "legend(\"topright\", c(\"Reference Sample to Generate PCs\", \"Study Sample Projected to Reference PC Space\"),\n");
      fprintf(fp, "col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = 19, cex = 1.2)\n");
   }else
      fprintf(fp, "plot(data$PC3, data$PC4, type=\"p\", xlab=\"PC3\", ylab=\"PC4\", main = \"Population Structure in %s\")\n",(const char*)prefix);
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   if (rpath)
       sprintf(command, "%s CMD BATCH %s", rpath, (const char*)scriptfile);
   else
       sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("--mds or --pca is done but R code %s failed.\n\n", (const char*)scriptfile);
   else if(error){
      printf("--ibdmap or --pca is done but R code %s failed.\n", (const char*)scriptfile);
      printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
   }else{
      sprintf(command, "ps2pdf %s_pcplot.ps", (const char*)prefix);
      system(command);
      printf("Population structure plot is generated in %s_pcplot.pdf\n", (const char*)prefix);
   }
}


int CheckRout(const char *scriptfile)
{
   String outfile=scriptfile;
   outfile.Add("out");
   int errorFlag = 1;
   char buff[256];   // define the buffer and allocate the length
   FILE *fp = fopen((const char*)outfile, "rb");
   if(fp != NULL){
      fseek(fp, -17, SEEK_END);// set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
      fread(buff, 16, 1, fp); // read the contents of the file starting from where fseek() positioned us
      buff[16] = '\0';        // close the string
      String lastline=buff;
      if(lastline.Find("Execution halted")==-1) errorFlag = 0;
      else{
         fseek(fp, -100, SEEK_END);  // set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
            fread(buff, 80, 1, fp); // read the contents of the file starting from where fseek() positioned us
         buff[80] = '\0';           // close the string
         lastline=buff;
         if(lastline.Find("Error in library")>-1) errorFlag = 3;
      }
      fclose(fp);                   // close the file
   }else errorFlag = 2;             // unable to open .Rout file
   return errorFlag;
}

void plotIBD1vsIBD2(const char *prefix, FILE *fp)
{
   fprintf(fp, "if(length(data$IBD1Seg)>0 & length(data$IBD2Seg)>0){\n");
   fprintf(fp, "d0 <- data$IBD2Seg>0.7\n");
   fprintf(fp, "d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)\n");
   fprintf(fp, "d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08\n");
   fprintf(fp, "d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)\n");
   fprintf(fp, "d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678\n");
   fprintf(fp, "d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839\n");
   fprintf(fp, "dU <- data$PropIBD>0 & data$PropIBD<=0.04419\n");
   fprintf(fp, "plot(data$IBD1Seg[dU], data$IBD2Seg[dU], type=\"p\", col = \"black\", cex.lab=1.2,\n");
   fprintf(fp, "xlim=c(min(data$IBD1Seg), max(data$IBD1Seg)),\n");
   fprintf(fp, "ylim=c(min(data$IBD2Seg), max(data$IBD2Seg)),\n");
   fprintf(fp, "main = \"IBD Segments In Inferred %s Relatives\",\n", (const char*)prefix);
   fprintf(fp, "xlab=expression(paste(\"Length Proportion of IBD1 Segments (\", pi[1], \")\",sep=\"\")),\n");
   fprintf(fp, "ylab=expression(paste(\"Length Proportion of IBD2 Segments (\", pi[2], \")\",sep=\"\")))\n");
   fprintf(fp, "points(data$IBD1Seg[d0], data$IBD2Seg[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.PO], data$IBD2Seg[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$IBD1Seg[d1.FS], data$IBD2Seg[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$IBD1Seg[d2], data$IBD2Seg[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$IBD1Seg[d3], data$IBD2Seg[d3], col=\"magenta\")\n");
   fprintf(fp, "points(data$IBD1Seg[d4], data$IBD2Seg[d4], col=\"gold\")\n");
   fprintf(fp, "abline(h = 0.08, col = \"green\", lty = 3, lwd = 2)\n");   // FS && MZ
   fprintf(fp, "abline(a = 0.96, b = -1, col = \"red\", lty = 3, lwd = 2)\n");   // PO
   fprintf(fp, "abline(a = 0.3535534, b = -0.5, col = \"green\", lty = 3, lwd = 2)\n");// FS
   fprintf(fp, "abline(a = 0.1767767, b = -0.5, col = \"blue\", lty = 3, lwd = 2)\n");// 2nd/3rd
   fprintf(fp, "abline(a = 0.08838835, b = -0.5, col = \"magenta\", lty = 3, lwd = 2)\n");// 3rd/4th
   fprintf(fp, "abline(a = 0.04419, b = -0.5, col = \"gold\", lty = 3, lwd = 2)\n");// 4th/UN
   fprintf(fp, "allcolors <- c(\"purple\", \"red\", \"green\", \"blue\", \"magenta\", \"gold\", \"black\")\n");
   fprintf(fp, "legend(\"topright\", c(\"Inferred MZ\", \"Inferred PO\", \"Inferred FS\", \"Inferred 2nd\", \"Inferred 3rd\", \"Inferred 4th\", \"Inferred UN\"),\n");
   fprintf(fp, "  col=allcolors, text.col = allcolors, pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
}


void MakeHashKG(StringIntHash & HashKG)
{
   const int POPCOUNT = 6;
   StringArray KG[POPCOUNT];
   String temp;
   temp = "HG01879 HG01880 HG01882 HG01883 HG01885 HG01956 HG01886 HG02014 HG01889 HG01890 HG01894 HG01912 HG01896 HG02013 HG01985 HG01914 HG01915 HG01958 HG01986 HG01988 HG01989 HG01990 HG02012 HG02009 HG02010 HG02051 HG02052 HG02053 HG02054 HG02095 HG02332 HG02107 HG02108 HG02111 HG02143 HG02144 HG02255 HG02256 HG02314 HG02315 HG02307 HG02308 HG02281 HG02282 HG02283 HG02309 HG02284 HG02317 HG02318 HG02322 HG02323 HG02325 HG02337 HG02334 HG02330 HG02339 HG02343 HG02419 HG02420 HG02427 HG02429 HG02433 HG02439 HG02442 HG02445 HG02489 HG02449 HG02450 HG02455 HG02470 HG02471 HG02476 HG02477 HG02479 HG02481 HG02484 HG02485 HG02511 HG02496 HG02497 HG02501 HG02502 HG02505 HG02508 HG02536 HG02537 HG02541 HG02545 HG02546 HG02549 HG02554 HG02555 HG02557 HG02558 HG02577 HG02580 NA19625 NA20274 NA19700 NA19701";
   temp += " NA19703 NA19704 NA19707 NA19711 NA19712 NA19818 NA19819 NA19834 NA19835 NA19900 NA19901 NA19904 NA19913 NA19908 NA19909 NA19914 NA19916 NA19917 NA19920 NA19921 NA19922 NA19923 NA19713 NA19982 NA19984 NA20126 NA20127 NA20276 NA20278 NA20282 NA20281 NA20287 NA20289 NA20291 NA20294 NA20296 NA20299 NA20298 NA20314 NA20317 NA20318 NA20321 NA20320 NA20332 NA20334 NA20355 NA20339 NA20340 NA20342 NA20346 NA20348 NA20351 NA20356 NA20357 NA20359 NA20362 NA20412 HG02922 HG02923 HG03499 HG03511 HG03514 HG03515 HG03517 HG03518 HG03520 HG03521 HG02938 HG02941 HG02943 HG02944 HG02946 HG02947 HG02952 HG02953 HG02968 HG02970 HG02971 HG02973 HG02974 HG02976 HG02977 HG02979 HG02981 HG03099 HG03100 HG03103 HG03105 HG03108 HG03109 HG03111 HG03112 HG03114 HG03115 HG03117 HG03118 HG03120 HG03121 HG03123 HG03124";
   temp += " HG03126 HG03127 HG03129 HG03130 HG03132 HG03133 HG03135 HG03136 HG03139 HG03157 HG03159 HG03160 HG03162 HG03163 HG03166 HG03168 HG03169 HG03172 HG03175 HG03189 HG03190 HG03193 HG03195 HG03196 HG03198 HG03199 HG03202 HG03265 HG03267 HG03268 HG03270 HG03271 HG03279 HG03280 HG03291 HG03294 HG03295 HG03297 HG03298 HG03300 HG03301 HG03303 HG03304 HG03311 HG03313 HG03342 HG03343 HG03351 HG03352 HG03354 HG03363 HG03366 HG03367 HG03369 HG03370 HG03372 HG03024 HG03025 HG03027 HG03028 HG03039 HG03040 HG03045 HG03046 HG03048 HG03049 HG03240 HG03241 HG03246 HG03247 HG03258 HG03259 HG03538 HG03539 HG02461 HG02462 HG02464 HG02465 HG02561 HG02562 HG02568 HG02570 HG02571 HG02573 HG02574 HG02582 HG02583 HG02585 HG02586 HG02588 HG02589 HG02594 HG02595 HG02610 HG02611 HG02613 HG02614 HG02620 HG02621 HG02623";
   temp += " HG02624 HG02628 HG02629 HG02634 HG02635 HG02642 HG02643 HG02645 HG02646 HG02666 HG02667 HG02675 HG02676 HG02678 HG02679 HG02702 HG02703 HG02715 HG02716 HG02721 HG02722 HG02756 HG02757 HG02759 HG02760 HG02763 HG02768 HG02769 HG02771 HG02772 HG02798 HG02799 HG02804 HG02805 HG02807 HG02808 HG02810 HG02811 HG02813 HG02814 HG02816 HG02817 HG02819 HG02820 HG02836 HG02837 HG02839 HG02840 HG02851 HG02852 HG02854 HG02855 HG02860 HG02861 HG02870 HG02878 HG02879 HG02881 HG02882 HG02884 HG02885 HG02887 HG02888 HG02890 HG02891 HG02895 HG02896 HG02982 HG02983 NA19331 NA19434 NA19445 NA19017 NA19019 NA19020 NA19023 NA19024 NA19025 NA19026 NA19027 NA19028 NA19030 NA19031 NA19035 NA19036 NA19037 NA19038 NA19041 NA19042 NA19043 NA19307 NA19308 NA19309 NA19310 NA19312 NA19314 NA19315 NA19316 NA19317 NA19318";
   temp += " NA19319 NA19320 NA19321 NA19323 NA19324 NA19327 NA19328 NA19332 NA19334 NA19338 NA19346 NA19347 NA19350 NA19351 NA19355 NA19360 NA19372 NA19374 NA19375 NA19376 NA19377 NA19378 NA19379 NA19380 NA19383 NA19384 NA19385 NA19390 NA19391 NA19393 NA19394 NA19395 NA19397 NA19399 NA19401 NA19403 NA19404 NA19428 NA19429 NA19430 NA19431 NA19435 NA19436 NA19437 NA19438 NA19439 NA19440 NA19443 NA19446 NA19448 NA19449 NA19451 NA19452 NA19454 NA19455 NA19456 NA19457 NA19461 NA19462 NA19463 NA19466 NA19467 NA19468 NA19471 NA19472 NA19473 NA19474 NA19475 HG03057 HG03060 HG03074 HG03077 HG03078 HG03084 HG03085 HG03086 HG03378 HG03397 HG03432 HG03433 HG03436 HG03437 HG03442 HG03445 HG03446 HG03457 HG03460 HG03461 HG03472 HG03547 HG03548 HG03556 HG03557 HG03558 HG03565 HG03567 HG03052 HG03054 HG03055 HG03058";
   temp += " HG03061 HG03063 HG03064 HG03066 HG03069 HG03072 HG03073 HG03079 HG03081 HG03082 HG03088 HG03091 HG03095 HG03096 HG03097 HG03209 HG03212 HG03224 HG03225 HG03376 HG03380 HG03382 HG03385 HG03388 HG03391 HG03394 HG03401 HG03410 HG03419 HG03428 HG03439 HG03449 HG03451 HG03452 HG03455 HG03458 HG03464 HG03470 HG03469 HG03473 HG03476 HG03478 HG03479 HG03484 HG03485 HG03559 HG03563 HG03571 HG03572 HG03575 HG03577 HG03578 HG03583 NA18486 NA18488 NA18489 NA18498 NA18499 NA18501 NA18502 NA18504 NA18505 NA19107 NA19108 NA18867 NA18868 NA18507 NA18508 NA18511 NA18510 NA18858 NA18516 NA18517 NA18519 NA18520 NA18522 NA18523 NA18870 NA18871 NA18853 NA18873 NA18874 NA18876 NA18877 NA18878 NA18879 NA18881 NA18856 NA18861 NA18907 NA18908 NA18864 NA18865 NA18909 NA18910 NA18912 NA18915 NA18916 NA18917 NA18923";
   temp += " NA18924 NA19197 NA19198 NA18933 NA18934 NA19184 NA19185 NA19092 NA19093 NA19095 NA19096 NA19102 NA19137 NA19138 NA19175 NA19200 NA19201 NA19171 NA19172 NA19204 NA19209 NA19210 NA19206 NA19207 NA19159 NA19160 NA19225 NA19222 NA19223 NA19116 NA19119 NA19121 NA19141 NA19152 NA19153 NA19149 NA19143 NA19144 NA19146 NA19147 NA19129 NA19113 NA19114 NA19256 NA19257 NA19117 NA19118 NA19130 NA19131 NA19098 NA19099 NA19213 NA19214 NA19189 NA19190 NA19235 NA19236 NA19238 NA19239 NA19247 NA19248";
   KG[0].AddTokens(temp, ' ');   // AFR
   temp = "HG01119 HG01121 HG01122 HG01112 HG01113 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01275 HG01277 HG01280 HG01281 HG01284 HG01341 HG01342 HG01344 HG01345 HG01348 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01362 HG01363 HG01365 HG01366 HG01369 HG01372 HG01374 HG01375 HG01377 HG01378 HG01383 HG01384 HG01389 HG01390 HG01431 HG01432 HG01435 HG01437 HG01438 HG01441 HG01440 HG01443 HG01444 HG01447 HG01455 HG01456 HG01459 HG01461 HG01462 HG01464 HG01465 HG01468 HG01474 HG01479 HG01485 HG01486 HG01488 HG01489 HG01491 HG01492 HG01494 HG01495 HG01497 HG01498 HG01550 HG01551 HG01556 NA19734 NA19735 NA19740 NA19741 NA19752 NA19764";
   temp += " NA19792 NA19648 NA19649 NA19669 NA19670 NA19676 NA19651 NA19652 NA19654 NA19655 NA19657 NA19658 NA19678 NA19679 NA19681 NA19682 NA19661 NA19684 NA19663 NA19664 NA19716 NA19717 NA19719 NA19720 NA19722 NA19723 NA19725 NA19726 NA19728 NA19729 NA19731 NA19732 NA19746 NA19747 NA19749 NA19750 NA19755 NA19756 NA19758 NA19759 NA19761 NA19762 NA19770 NA19771 NA19785 NA19786 NA19773 NA19774 NA19776 NA19777 NA19779 NA19780 NA19782 NA19783 NA19788 NA19789 NA19794 NA19795 HG01565 HG01566 HG01571 HG01572 HG01577 HG01578 HG01917 HG01918 HG01892 HG01893 HG01920 HG01921 HG01923 HG01924 HG01926 HG01927 HG01932 HG01933 HG01935 HG01936 HG01938 HG01939 HG01941 HG01942 HG01944 HG01945 HG01947 HG01948 HG01950 HG01951 HG01953 HG01954 HG01967 HG01968 HG01970 HG01971 HG01961 HG01965 HG01973 HG01974 HG01976 HG01977";
   temp += " HG01979 HG01980 HG01982 HG01991 HG01992 HG01997 HG02008 HG02002 HG02003 HG02006 HG02089 HG02090 HG02102 HG02104 HG02105 HG02146 HG02147 HG02150 HG02252 HG02253 HG02259 HG02260 HG02262 HG02265 HG02266 HG02271 HG02272 HG02274 HG02275 HG02277 HG02278 HG02285 HG02286 HG02291 HG02292 HG02298 HG02299 HG02301 HG02304 HG02312 HG02345 HG02348 HG02425 HG00551 HG01075 HG00553 HG00554 HG00637 HG00638 HG00640 HG00641 HG00731 HG00732 HG00734 HG01047 HG00736 HG00737 HG01048 HG01049 HG00739 HG00740 HG00742 HG00743 HG01051 HG01052 HG01054 HG01055 HG01058 HG01060 HG01061 HG01063 HG01064 HG01066 HG01067 HG01069 HG01070 HG01072 HG01073 HG01101 HG01102 HG01077 HG01079 HG01080 HG01241 HG01242 HG01170 HG01171 HG01104 HG01105 HG01085 HG01086 HG01088 HG01089 HG01082 HG01083 HG01094 HG01095 HG01092 HG01097 HG01098";
   temp += " HG01247 HG01248 HG01107 HG01108 HG01110 HG01111 HG01161 HG01162 HG01173 HG01174 HG01164 HG01187 HG01188 HG01167 HG01168 HG01190 HG01191 HG01176 HG01177 HG01182 HG01183 HG01197 HG01198 HG01200 HG01204 HG01205 HG01286 HG01302 HG01303 HG01305 HG01308 HG01311 HG01312 HG01323 HG01325 HG01326 HG01392 HG01393 HG01395 HG01396 HG01398 HG01402 HG01403 HG01405 HG01412 HG01413 HG01414";
   KG[1].AddTokens(temp, ' ');   // AMR
   temp = "HG00867 HG02371 HG00759 HG00766 HG00844 HG00851 HG00864 HG00879 HG00881 HG00956 HG00978 HG00982 HG01028 HG01029 HG01031 HG01046 HG01794 HG01795 HG01796 HG01797 HG01798 HG01799 HG01800 HG01801 HG01802 HG01804 HG01805 HG01806 HG01807 HG01808 HG01809 HG01810 HG01811 HG01812 HG01813 HG01815 HG01816 HG01817 HG02151 HG02152 HG02153 HG02154 HG02155 HG02156 HG02164 HG02165 HG02166 HG02178 HG02179 HG02180 HG02181 HG02182 HG02184 HG02185 HG02186 HG02187 HG02188 HG02190 HG02250 HG02351 HG02353 HG02355 HG02356 HG02360 HG02364 HG02367 HG02373 HG02374 HG02375 HG02379 HG02380 HG02382 HG02383 HG02384 HG02385 HG02386 HG02389 HG02390 HG02391 HG02392 HG02394 HG02395 HG02396 HG02397 HG02398 HG02399 HG02401 HG02402 HG02406 HG02407 HG02408 HG02409 HG02410 NA18525 NA18526 NA18528 NA18530 NA18531 NA18532 NA18533";
   temp += " NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18570 NA18571 NA18572 NA18573 NA18574 NA18577 NA18579 NA18582 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 NA18599 NA18602 NA18603 NA18605 NA18606 NA18608 NA18609 NA18610 NA18611 NA18612 NA18613 NA18614 NA18615 NA18616 NA18617 NA18618 NA18619 NA18620 NA18621 NA18622 NA18623 NA18624 NA18625 NA18626 NA18627 NA18628 NA18629 NA18630 NA18631 NA18632 NA18633 NA18634 NA18635 NA18636 NA18637 NA18638 NA18639 NA18640 NA18641 NA18642 NA18643 NA18644 NA18645 NA18646 NA18647 NA18648 NA18740 NA18745 NA18747 NA18748 NA18749 NA18757 HG00403 HG00404 HG00406 HG00407";
   temp += " HG00409 HG00410 HG00419 HG00421 HG00422 HG00428 HG00436 HG00437 HG00442 HG00443 HG00445 HG00446 HG00448 HG00449 HG00451 HG00452 HG00457 HG00458 HG00463 HG00464 HG00472 HG00473 HG00475 HG00476 HG00478 HG00479 HG00500 HG00513 HG00524 HG00525 HG00530 HG00531 HG00533 HG00534 HG00536 HG00537 HG00542 HG00543 HG00556 HG00557 HG00559 HG00560 HG00565 HG00566 HG00580 HG00581 HG00583 HG00584 HG00589 HG00590 HG00592 HG00593 HG00595 HG00596 HG00598 HG00599 HG00607 HG00608 HG00610 HG00611 HG00613 HG00614 HG00619 HG00620 HG00622 HG00623 HG00625 HG00626 HG00628 HG00629 HG00631 HG00632 HG00634 HG00650 HG00651 HG00653 HG00654 HG00662 HG00663 HG00671 HG00672 HG00674 HG00675 HG00683 HG00684 HG00689 HG00690 HG00692 HG00693 HG00698 HG00699 HG00656 HG00657 HG00701 HG00704 HG00705 HG00707 HG00708 HG00717 HG00728";
   temp += " HG00729 NA18939 NA18940 NA18941 NA18942 NA18943 NA18944 NA18945 NA18946 NA18947 NA18948 NA18949 NA18950 NA18951 NA18952 NA18953 NA18954 NA18956 NA18957 NA18959 NA18960 NA18961 NA18962 NA18963 NA18964 NA18965 NA18966 NA18967 NA18968 NA18969 NA18970 NA18971 NA18972 NA18973 NA18974 NA18975 NA18976 NA18977 NA18978 NA18979 NA18980 NA18981 NA18982 NA18983 NA18984 NA18985 NA18986 NA18987 NA18988 NA18989 NA18990 NA18991 NA18992 NA18993 NA18994 NA18995 NA18997 NA18998 NA18999 NA19000 NA19001 NA19002 NA19003 NA19004 NA19005 NA19006 NA19007 NA19009 NA19010 NA19011 NA19012 NA19054 NA19055 NA19056 NA19057 NA19058 NA19059 NA19060 NA19062 NA19063 NA19064 NA19065 NA19066 NA19067 NA19068 NA19070 NA19072 NA19074 NA19075 NA19076 NA19077 NA19078 NA19079 NA19080 NA19081 NA19082 NA19083 NA19084 NA19085 NA19086";
   temp += " NA19087 NA19088 NA19089 NA19090 NA19091 HG01595 HG01596 HG01597 HG01598 HG01599 HG01600 HG01840 HG01841 HG01842 HG01843 HG01844 HG01845 HG01846 HG01847 HG01848 HG01849 HG01850 HG01851 HG01852 HG01853 HG01855 HG01857 HG01858 HG01859 HG01860 HG01861 HG01862 HG01863 HG01864 HG01865 HG01866 HG01867 HG01868 HG01869 HG01870 HG01871 HG01872 HG01873 HG01874 HG01878 HG02016 HG02017 HG02019 HG02020 HG02023 HG02025 HG02026 HG02028 HG02029 HG02031 HG02032 HG02035 HG02040 HG02047 HG02048 HG02049 HG02050 HG02057 HG02058 HG02060 HG02061 HG02064 HG02067 HG02069 HG02070 HG02072 HG02073 HG02075 HG02076 HG02078 HG02079 HG02081 HG02082 HG02084 HG02085 HG02086 HG02087 HG02088 HG02113 HG02116 HG02121 HG02122 HG02127 HG02128 HG02130 HG02131 HG02133 HG02134 HG02136 HG02137 HG02138 HG02139 HG02140 HG02141 HG02142 HG02512 HG02513 HG02521 HG02522";
   KG[2].AddTokens(temp, ' ');   // EAS
   temp = "NA06984 NA06989 NA12347 NA12348 NA06986 NA07037 NA07051 NA12340 NA12341 NA12342 NA10847 NA12144 NA06994 NA07000 NA07056 NA06985 NA07048 NA10851 NA12058 NA07347 NA07357 NA12043 NA12044 NA12045 NA12046 NA11881 NA11840 NA11843 NA11829 NA11830 NA11831 NA11832 NA12383 NA12489 NA12546 NA12399 NA12400 NA12413 NA12414 NA12716 NA12717 NA12718 NA11992 NA11994 NA11995 NA12234 NA11892 NA11893 NA11894 NA12154 NA12155 NA12156 NA12249 NA12272 NA12273 NA12275 NA12003 NA12004 NA12005 NA12006 NA12286 NA12287 NA12282 NA12283 NA11918 NA11919 NA11920 NA11930 NA11931 NA11932 NA11933 NA12748 NA12749 NA12750 NA12751 NA12760 NA12761 NA12762 NA12763 NA12775 NA12776 NA12777 NA12778 NA12812 NA12813 NA12814 NA12815 NA12827 NA12828 NA12829 NA12830 NA12842 NA12843 NA12872 NA12873 NA12874 NA12878 NA12889 NA12890 HG00171";
   temp += " HG00173 HG00174 HG00176 HG00177 HG00178 HG00179 HG00180 HG00181 HG00182 HG00183 HG00185 HG00186 HG00187 HG00188 HG00189 HG00190 HG00266 HG00267 HG00268 HG00269 HG00271 HG00272 HG00273 HG00274 HG00275 HG00276 HG00277 HG00278 HG00280 HG00281 HG00282 HG00284 HG00285 HG00288 HG00290 HG00304 HG00306 HG00308 HG00309 HG00310 HG00311 HG00313 HG00315 HG00318 HG00319 HG00320 HG00321 HG00323 HG00324 HG00325 HG00326 HG00327 HG00328 HG00329 HG00330 HG00331 HG00332 HG00334 HG00335 HG00336 HG00337 HG00338 HG00339 HG00341 HG00342 HG00343 HG00344 HG00345 HG00346 HG00349 HG00350 HG00351 HG00353 HG00355 HG00356 HG00357 HG00358 HG00360 HG00361 HG00362 HG00364 HG00365 HG00366 HG00367 HG00368 HG00369 HG00371 HG00372 HG00373 HG00375 HG00376 HG00378 HG00379 HG00380 HG00381 HG00382 HG00383 HG00384 HG00155 HG00146";
   temp += " HG00158 HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00105 HG00106 HG00107 HG00108 HG00109 HG00110 HG00111 HG00112 HG00113 HG00114 HG00115 HG00116 HG00117 HG00118 HG00119 HG00120 HG00121 HG00122 HG00123 HG00125 HG00126 HG00127 HG00128 HG00129 HG00130 HG00131 HG00132 HG00133 HG00136 HG00137 HG00138 HG00139 HG00140 HG00141 HG00142 HG00143 HG00145 HG00148 HG00149 HG00150 HG00151 HG00154 HG00157 HG00159 HG00160 HG00231 HG00232 HG00233 HG00234 HG00235 HG00236 HG00237 HG00238 HG00239 HG00240 HG00242 HG00243 HG00244 HG00245 HG00246 HG00250 HG00251 HG00252 HG00253 HG00254 HG00255 HG00256 HG00257 HG00258 HG00259 HG00260 HG00261 HG00262 HG00263 HG00264 HG00265 HG01334 HG01789 HG01790 HG01791 HG02215 HG01500 HG01501 HG01503 HG01504 HG01506 HG01507 HG01509 HG01510 HG01512 HG01513 HG01515";
   temp += " HG01516 HG01518 HG01519 HG01521 HG01522 HG01524 HG01525 HG01527 HG01528 HG01530 HG01531 HG01536 HG01537 HG01631 HG01632 HG01628 HG01630 HG01625 HG01626 HG01623 HG01624 HG01619 HG01620 HG01617 HG01618 HG01613 HG01615 HG01610 HG01612 HG01607 HG01608 HG01605 HG01606 HG01602 HG01603 HG01668 HG01669 HG01670 HG01672 HG01673 HG01675 HG01676 HG01678 HG01679 HG01680 HG01682 HG01684 HG01685 HG01686 HG01694 HG01695 HG01697 HG01699 HG01700 HG01702 HG01704 HG01705 HG01707 HG01708 HG01709 HG01710 HG01746 HG01747 HG01756 HG01757 HG01761 HG01762 HG01765 HG01766 HG01767 HG01768 HG01770 HG01771 HG01773 HG01775 HG01776 HG01777 HG01779 HG01781 HG01783 HG01784 HG01785 HG01786 HG02219 HG02220 HG02221 HG02223 HG02224 HG02230 HG02231 HG02232 HG02233 HG02235 HG02236 HG02238 HG02239 NA20502 NA20503 NA20504 NA20505";
   temp += " NA20506 NA20507 NA20508 NA20509 NA20510 NA20511 NA20512 NA20513 NA20514 NA20515 NA20516 NA20517 NA20518 NA20519 NA20520 NA20521 NA20522 NA20524 NA20525 NA20527 NA20528 NA20529 NA20530 NA20531 NA20532 NA20533 NA20534 NA20535 NA20536 NA20538 NA20539 NA20540 NA20541 NA20542 NA20543 NA20544 NA20581 NA20582 NA20585 NA20586 NA20587 NA20588 NA20589 NA20752 NA20753 NA20754 NA20755 NA20756 NA20757 NA20758 NA20759 NA20760 NA20761 NA20762 NA20763 NA20764 NA20765 NA20766 NA20767 NA20768 NA20769 NA20770 NA20771 NA20772 NA20773 NA20774 NA20775 NA20778 NA20783 NA20785 NA20786 NA20787 NA20790 NA20792 NA20795 NA20796 NA20797 NA20798 NA20799 NA20800 NA20801 NA20802 NA20803 NA20804 NA20805 NA20806 NA20807 NA20808 NA20809 NA20810 NA20811 NA20812 NA20813 NA20814 NA20815 NA20818 NA20819 NA20821 NA20822 NA20826 NA20827 NA20828 NA20832";
   KG[3].AddTokens(temp, ' ');   // EUR
   temp = "HG03006 HG03007 HG03589 HG03600 HG03603 HG03604 HG03607 HG03611 HG03615 HG03616 HG03793 HG03796 HG03800 HG03802 HG03803 HG03805 HG03808 HG03809 HG03814 HG03815 HG03817 HG03821 HG03823 HG03824 HG03826 HG03829 HG03830 HG03832 HG03833 HG03012 HG03905 HG03907 HG03908 HG03910 HG03911 HG03913 HG03914 HG03917 HG03919 HG03920 HG03922 HG03925 HG03926 HG03928 HG03931 HG03585 HG03934 HG03937 HG03940 HG03941 HG04131 HG04134 HG04140 HG04141 HG04144 HG04146 HG04152 HG04153 HG04155 HG04156 HG04158 HG04159 HG04161 HG04162 HG04164 HG04171 HG04173 HG04176 HG04177 HG03593 HG04180 HG04182 HG04183 HG04185 HG04186 HG04188 HG04189 HG04194 HG04195 HG03594 HG03595 HG03598 HG03009 HG03812 HG03902 HG03916 NA20868 NA20886 NA20910 NA20845 NA20846 NA20847 NA20849 NA20850 NA20851 NA20852 NA20853 NA20854 NA20856 NA20858";
   temp += " NA20859 NA20861 NA20862 NA20863 NA20864 NA20866 NA20867 NA20869 NA20870 NA20872 NA20874 NA20875 NA20876 NA20877 NA20878 NA20881 NA20882 NA20884 NA20885 NA20887 NA20888 NA20889 NA20890 NA20891 NA20892 NA20894 NA20895 NA20896 NA20897 NA20899 NA20900 NA20901 NA20902 NA20903 NA20904 NA20905 NA20906 NA20908 NA20911 NA21086 NA21087 NA21088 NA21089 NA21090 NA21091 NA21092 NA21093 NA21094 NA21095 NA21097 NA21098 NA21099 NA21100 NA21101 NA21102 NA21103 NA21104 NA21105 NA21106 NA21107 NA21108 NA21109 NA21110 NA21111 NA21112 NA21113 NA21114 NA21115 NA21116 NA21117 NA21118 NA21119 NA21120 NA21122 NA21123 NA21124 NA21125 NA21126 NA21127 NA21128 NA21129 NA21130 NA21133 NA21135 NA21137 NA21141 NA21142 NA21143 NA21144 HG03871 HG04206 HG04239 HG03713 HG03722 HG03727 HG03772 HG03773 HG03874 HG04018 HG03714";
   temp += " HG03716 HG03717 HG03718 HG03730 HG03729 HG03731 HG03742 HG03770 HG03771 HG03781 HG03779 HG03774 HG03775 HG03777 HG03784 HG03786 HG03785 HG03790 HG03792 HG03787 HG03788 HG03789 HG03782 HG03866 HG03873 HG03861 HG03863 HG03864 HG03862 HG03869 HG03872 HG03870 HG03867 HG03780 HG03882 HG03963 HG03960 HG03968 HG03969 HG03977 HG03978 HG03971 HG03974 HG03973 HG03976 HG03720 HG04001 HG04002 HG04019 HG04023 HG04025 HG04014 HG04020 HG04015 HG04017 HG04026 HG04022 HG04054 HG04056 HG04061 HG04062 HG04060 HG04063 HG04059 HG04076 HG04080 HG04090 HG04093 HG04094 HG04096 HG04098 HG03778 HG04118 HG04209 HG04198 HG04200 HG04202 HG04211 HG04212 HG04222 HG04214 HG04216 HG04219 HG04225 HG04235 HG04238 HG03875 HG03965 HG04070 HG03868 HG03967 HG03022 HG01583 HG01586 HG01589 HG01593 HG02490 HG02491 HG02493 HG02494";
   temp += " HG02597 HG02600 HG02601 HG02603 HG02604 HG02648 HG02649 HG02651 HG02652 HG02654 HG02655 HG02657 HG02658 HG02660 HG02661 HG02681 HG02682 HG02684 HG02685 HG02687 HG02688 HG02690 HG02691 HG02694 HG02696 HG02697 HG02699 HG02700 HG02724 HG02725 HG02727 HG02728 HG02731 HG02733 HG02734 HG02736 HG02737 HG02774 HG02775 HG02778 HG02780 HG02783 HG02784 HG02786 HG02787 HG02789 HG02790 HG02792 HG02793 HG03015 HG03016 HG03018 HG03019 HG03021 HG03228 HG03229 HG03234 HG03235 HG03237 HG03238 HG03488 HG03490 HG03491 HG03619 HG03624 HG03625 HG03629 HG03631 HG03634 HG03636 HG03640 HG03649 HG03652 HG03653 HG03660 HG03663 HG03667 HG03668 HG03702 HG03703 HG03705 HG03706 HG03708 HG03709 HG03762 HG03765 HG03767 HG03646 HG03645 HG03643 HG03644 HG03642 HG03679 HG03999 HG03672 HG03680 HG03673 HG03681 HG03684 HG03685";
   temp += " HG03951 HG03686 HG03687 HG03885 HG03886 HG03689 HG03690 HG03691 HG03692 HG03693 HG03995 HG03694 HG03695 HG03696 HG03697 HG03698 HG03711 HG03740 HG03741 HG03884 HG03733 HG03736 HG03738 HG03743 HG03757 HG03744 HG03745 HG03746 HG03888 HG03750 HG03754 HG03755 HG03756 HG03752 HG03753 HG03760 HG03836 HG03844 HG04033 HG03838 HG04006 HG03837 HG03846 HG03849 HG03850 HG03848 HG03851 HG03858 HG03950 HG03857 HG03856 HG03854 HG04035 HG03887 HG03986 HG03895 HG03890 HG03894 HG03896 HG03898 HG03899 HG03897 HG03900 HG03943 HG03944 HG03945 HG03947 HG03955 HG03953 HG03949 HG03990 HG03991 HG03985 HG03989 HG03998 HG04003 HG04029 HG04038 HG04042 HG04039 HG04047 HG04106 HG04107 HG04210 HG04075 HG04099 HG04100 HG04227 HG04229";
   KG[4].AddTokens(temp, ' ');   // SAS
   temp = "HG00542 HG00739 HG01075 HG01108 HG01241 HG01275 HG01438 HG02479 HG03754 HG03899 HG03998 NA19042 NA19334 NA19679 NA19913 NA20274 NA20314 NA20318 NA20321 NA20355 NA20362 NA20900 NA21135";
   KG[5].AddTokens(temp, ' ');   // Excluded
   HashKG.Clear();
   for(int i = 0; i < POPCOUNT; i ++)
      for(int j = 0; j < KG[i].Length(); j++)
         HashKG.SetInteger(KG[i][j], i+1);
}

void MakeHashEUR(StringIntHash & HashEUR)
{
const int POPCOUNT=4;
StringArray EUR[POPCOUNT];
String temp = "NA06984 NA06989 NA12347 NA12348 NA06986 NA07037 NA07051 NA12340 NA12341 NA12342 NA10847 NA12144 NA06994 NA07000 NA07056 NA06985 NA07048 NA10851 NA12058 NA07347 NA07357 NA12043 NA12044 NA12045 NA12046 NA11881 NA11840 NA11843 NA11829 NA11830 NA11831 NA11832 NA12383 NA12489 NA12546 NA12399 NA12400 NA12413 NA12414 NA12716 NA12717 NA12718 NA11992 NA11994 NA11995 NA12234 NA11892 NA11893 NA11894 NA12154 NA12155 NA12156 NA12249 NA12272 NA12273 NA12275 NA12003 NA12004 NA12005 NA12006 NA12286 NA12287 NA12282 NA12283 NA11918 NA11919 NA11920 NA11930 NA11931 NA11932 NA11933 NA12748 NA12749 NA12750 NA12751 NA12760 NA12761 NA12762 NA12763 NA12775 NA12776 NA12777 NA12778 NA12812 NA12813 NA12814 NA12815 NA12827 NA12828 NA12829 NA12830 NA12842 NA12843 NA12872 NA12873 NA12874 NA12878 NA12889 NA12890 HG00155";
temp += " HG00146 HG00158 HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00105 HG00106 HG00107 HG00108 HG00109 HG00110 HG00111 HG00112 HG00113 HG00114 HG00115 HG00116 HG00117 HG00118 HG00119 HG00120 HG00121 HG00122 HG00123 HG00125 HG00126 HG00127 HG00128 HG00129 HG00130 HG00131 HG00132 HG00133 HG00136 HG00137 HG00138 HG00139 HG00140 HG00141 HG00142 HG00143 HG00145 HG00148 HG00149 HG00150 HG00151 HG00154 HG00157 HG00159 HG00160 HG00231 HG00232 HG00233 HG00234 HG00235 HG00236 HG00237 HG00238 HG00239 HG00240 HG00242 HG00243 HG00244 HG00245 HG00246 HG00250 HG00251 HG00252 HG00253 HG00254 HG00255 HG00256 HG00257 HG00258 HG00259 HG00260 HG00261 HG00262 HG00263 HG00264 HG00265 HG01334 HG01789 HG01790 HG01791 HG02215";
EUR[0].AddTokens(temp, ' ');   // NEUR
temp = "HG01500 HG01501 HG01503 HG01504 HG01506 HG01507 HG01509 HG01510 HG01512 HG01513 HG01515 HG01516 HG01518 HG01519 HG01521 HG01522 HG01524 HG01525 HG01527 HG01528 HG01530 HG01531 HG01536 HG01537 HG01631 HG01632 HG01628 HG01630 HG01625 HG01626 HG01623 HG01624 HG01619 HG01620 HG01617 HG01618 HG01613 HG01615 HG01610 HG01612 HG01607 HG01608 HG01605 HG01606 HG01602 HG01603 HG01668 HG01669 HG01670 HG01672 HG01673 HG01675 HG01676 HG01678 HG01679 HG01680 HG01682 HG01684 HG01685 HG01686 HG01694 HG01695 HG01697 HG01699 HG01700 HG01702 HG01704 HG01705 HG01707 HG01708 HG01709 HG01710 HG01746 HG01747 HG01756 HG01757 HG01761 HG01762 HG01765 HG01766 HG01767 HG01768 HG01770 HG01771 HG01773 HG01775 HG01776 HG01777 HG01779 HG01781 HG01783 HG01784 HG01785 HG01786 HG02219 HG02220 HG02221 HG02223 HG02224 HG02230";
temp += " HG02231 HG02232 HG02233 HG02235 HG02236 HG02238 HG02239 NA20502 NA20503 NA20504 NA20505 NA20506 NA20507 NA20508 NA20509 NA20510 NA20511 NA20512 NA20513 NA20514 NA20515 NA20516 NA20517 NA20518 NA20519 NA20520 NA20521 NA20522 NA20524 NA20525 NA20527 NA20528 NA20529 NA20530 NA20531 NA20532 NA20533 NA20534 NA20535 NA20536 NA20538 NA20539 NA20540 NA20541 NA20542 NA20543 NA20544 NA20581 NA20582 NA20585 NA20586 NA20587 NA20588 NA20589 NA20752 NA20753 NA20754 NA20755 NA20756 NA20757 NA20758 NA20759 NA20760 NA20761 NA20762 NA20763 NA20764 NA20765 NA20766 NA20767 NA20768 NA20769 NA20770 NA20771 NA20772 NA20773 NA20774 NA20775 NA20778 NA20783 NA20785 NA20786 NA20787 NA20790 NA20792 NA20795 NA20796 NA20797 NA20798 NA20799 NA20800 NA20801 NA20802 NA20803 NA20804 NA20805 NA20806 NA20807 NA20808 NA20809 NA20810 NA20811 NA20812 NA20813 NA20814 NA20815 NA20818 NA20819 NA20821 NA20822 NA20826 NA20827 NA20828 NA20832";
EUR[1].AddTokens(temp, ' ');   // SEUR
temp = "HG00171 HG00173 HG00174 HG00176 HG00177 HG00178 HG00179 HG00180 HG00181 HG00182 HG00183 HG00185 HG00186 HG00187 HG00188 HG00189 HG00190 HG00266 HG00267 HG00268 HG00269 HG00271 HG00272 HG00273 HG00274 HG00275 HG00276 HG00277 HG00278 HG00280 HG00281 HG00282 HG00284 HG00285 HG00288 HG00290 HG00304 HG00306 HG00308 HG00309 HG00310 HG00311 HG00313 HG00315 HG00318 HG00319 HG00320 HG00321 HG00323 HG00324 HG00325 HG00326 HG00327 HG00328 HG00329 HG00330 HG00331 HG00332 HG00334 HG00335 HG00336 HG00337 HG00338 HG00339 HG00341 HG00342 HG00343 HG00344 HG00345 HG00346 HG00349 HG00351 HG00353 HG00355 HG00356 HG00357 HG00358 HG00360 HG00361 HG00362 HG00364 HG00365 HG00366 HG00367 HG00368 HG00369 HG00371 HG00372 HG00373 HG00375 HG00376 HG00378 HG00379 HG00380 HG00381 HG00382 HG00383 HG00384";
EUR[2].AddTokens(temp, ' ');   // FIN
HashEUR.Clear();
temp = "HG00350";
EUR[3].AddTokens(temp, ' ');   // Excluded
HashEUR.Clear();
for(int i = 0; i < POPCOUNT; i ++)
   for(int j = 0; j < EUR[i].Length(); j++)
      HashEUR.SetInteger(EUR[i][j], i+1);
}

/*
void plotHetConcvsIBD2(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd2plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "#%s for KING, by Wei-Min Chen\n", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd2plot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%s.seg2\", header=T)\n", (const char*)prefix);
   fprintf(fp, "valid <- data$Pr_IBD2>0.005\n");
   fprintf(fp, "if(sum(valid)>0){\n");
   fprintf(fp, "plot(data$Pr_IBD2[valid], data$HetConc[valid], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "main = \"Relationships In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of IBD2 Segments\", ylab = \"Heterozygote Concordance Rate\")\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("R code %s cannot be run.\n", (const char*)scriptfile);
   else if(error)
      printf("Errors found in R code %s. Please check %sout for details.\n\n", (const char*)scriptfile, (const char*)scriptfile);
   else{
      sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
      system(command);
      printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
   }
}

void plotIBD2(const char *prefix)
{
   String scriptfile=prefix;
   scriptfile.Add("_ibd2plot.R");
   FILE *fp = fopen(scriptfile, "wt");
   if(fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "postscript(\"%s_ibd2plot.ps\", paper=\"letter\", horizontal=T)\n",
      (const char*)prefix);
   fprintf(fp, "data <- read.table(file=\"%s.ibs\", header=T)\n", (const char*)prefix);
   fprintf(fp, "if(dim(data)[1]>0){\n");
   fprintf(fp, "d0 <- data$Phi==0.5\n");
   fprintf(fp, "d1.PO <- data$Phi==0.25 & data$Z0==0\n");
   fprintf(fp, "d1.FS <- data$Phi==0.25 & data$Z0>0\n");
   fprintf(fp, "d2 <- data$Phi==0.125\n");
   fprintf(fp, "dU <- data$Phi==0\n");
   fprintf(fp, "dO <- !d0 & !d1.PO & !d1.FS & !d2 & !dU\n");
   fprintf(fp, "plot(data$Pr_IBD2[dU], data$Kinship[dU], type=\"p\",\n");
   fprintf(fp, "col = \"black\", cex.lab=1.3,\n");
   fprintf(fp, "xlim = c(min(data$Pr_IBD2), max(data$Pr_IBD2)),\n");
   fprintf(fp, "ylim = c(min(data$Kinship), max(data$Kinship)),\n");
   fprintf(fp, "main = \"Relationships In %s Families\",\n", (const char*)prefix);
   fprintf(fp, "xlab=\"Proportion of IBD2 Segments\", ylab = \"Estimated Kinship Coefficient\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d0], data$Kinship[d0], col=\"purple\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.PO], data$Kinship[d1.PO], col=\"red\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.FS], data$Kinship[d1.FS], col=\"green\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d2], data$Kinship[d2], col=\"blue\")\n");
   fprintf(fp, "points(data$Pr_IBD2[dO], data$Kinship[dO], col=\"gold\")\n");
   fprintf(fp, "points(data$Pr_IBD2[dU & data$Kinship>0.088], data$Kinship[dU & data$Kinship>0.088], col=\"black\")\n");
   fprintf(fp, "points(data$Pr_IBD2[d1.FS & data$Pr_IBD2<0.05], data$Kinship[d1.FS & data$Pr_IBD2<0.05], col=\"green\")\n");
   fprintf(fp, "abline(h = 0.3536, col = \"red\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.1768, col = \"green\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0884, col = \"blue\", lty = 3)\n");
   fprintf(fp, "abline(h = 0.0442, col = \"black\", lty = 3)\n");
   fprintf(fp, "legend(\"bottomright\", c(\"MZ Twin\", \"Parent-Offspring\", \"Full Siblings\", \"2nd-Degree\", \"More Distant\", \"Unrelated\"),\n");
   fprintf(fp, "col=c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"),\n");
   fprintf(fp, "text.col = c(\"purple\", \"red\", \"green\", \"blue\", \"gold\", \"black\"), pch = 19, cex = 1.2)\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fclose(fp);
   char command[256];
   sprintf(command, "R CMD BATCH %s", (const char*)scriptfile);
   system(command);
   int error = CheckRout(scriptfile);
   if(error == 2)
      printf("R code %s cannot be run.\n", (const char*)scriptfile);
   else if(error)
      printf("Errors found in R code %s. Please check %sout for details.\n\n", (const char*)scriptfile, (const char*)scriptfile);
   else{
      sprintf(command, "ps2pdf %s_ibd2plot.ps", (const char*)prefix);
      system(command);
      printf("  Relationship plot is generated in %s_ibd2plot.pdf\n", (const char*)prefix);
   }
}
*/

