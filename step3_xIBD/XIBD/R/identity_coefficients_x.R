# Internal Function
#
# Calculates the condensed identity coefficients for the X chromosome and
# then calculate IBD coefficients. The original formulation of the mathematics
# was done by David Wakeham.

partitions <- function(pedigree,n,i,j){
  sex_i <- pedigree[pedigree[,"iid"] == i,"sex"]
  sex_j <- pedigree[pedigree[,"iid"] == j,"sex"]
  Part <- 1
  Coef <- 1
  if(sex_i==1&sex_j==1){
    if(n==1){
      Block <- 1; Ind <- c(i,j); Allele <- 1
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==2){
      Block <- c(1,2); Ind <- c(i,j); Allele <- 1
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
  }
  if(sex_i==1&sex_j==2){
    if(n==1){
      Block <- 1; Ind <- c(i,j,j); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==2){
      Block <- c(1,1,2); Ind <- c(j,j,i); Allele <- c(1,2,1)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==3){
      Block <- c(1,1,2); Ind <- c(i,j,j); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==4){
      Block <- c(1,2,3); Ind <- c(i,j,j); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
  }
  if(sex_i==2&sex_j==1){
    if(n==1){
      Block <- 1; Ind <- c(j,i,i); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==2){
      Block <- c(1,1,2); Ind <- c(i,i,j); Allele <- c(1,2,1)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==3){
      Block <- c(1,1,2); Ind <- c(j,i,i); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==4){
      Block <- c(1,2,3); Ind <- c(j,i,i); Allele <- c(1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
  }
  if(sex_i==2&sex_j==2){
    if(n==1){
      Block <- 1; Ind <- c(i,i,j,j); Allele <- c(1,2,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==2){
      Block <- c(1,1,2,2); Ind <- c(i,i,j,j); Allele <- c(1,2,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==3){
      Block <- c(1,1,1,2); Ind <- c(i,i,j,j); Allele <- c(1,2,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==4){
      Block <- c(1,1,2,3); Ind <- c(i,i,j,j); Allele <- c(1,2,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==5){
      Block <- c(1,1,1,2); Ind <- c(i,j,j,i); Allele <- c(1,1,2,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==6){
      Block <- c(1,2,3,3); Ind <- c(i,i,j,j); Allele <- c(1,2,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==7){
      Block <- c(1,1,2,2); Ind <- c(i,j,i,j); Allele <- c(1,1,1,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==8){
      Block <- c(1,1,2,3); Ind <- c(i,j,i,j); Allele <- c(1,1,2,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
    if(n==9){
      Block <- c(1,2,3,4); Ind <- c(i,j,i,j); Allele <- c(1,1,2,2)
      partition.0 <- cbind(Part,Block,Ind,Allele,Coef)
    }
  }
  return(partition.0)
}


# selects the individual to perform the recurrence rules on
individual_i <- function(pedigree, partition){
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }

  founders <- pedigree[pedigree[,"pid"]==0&pedigree[,"mid"]==0,]

  for(found in founders[,"iid"]){
    partition <- matrix(partition[!partition[,"Ind"]==found,],ncol=5,dimnames=dimnames.partition)
  }
  return(max(partition[,"Ind"]))
}


recurrence_selection <- function(pedigree, partition){
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }

  recurrence <- NULL
  for(part in partition[!duplicated(partition[,"Part"]),"Part"]){
    partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
    ind.i <- individual_i(pedigree, partition.subset)
    blocks.w.i <- partition.subset[partition.subset[,"Ind"]==ind.i,"Block"][!duplicated(partition.subset[partition.subset[,"Ind"]==ind.i,"Block"])]
    if(length(blocks.w.i)==1){R <- 1 } else { R <- 2}
    P <- part
    recurrence <- rbind(recurrence, cbind(P,R))
  }
  return(recurrence)
}



recurrence_rule_1 <- function(pedigree, partition, recurrence){
  # Step 1. Check that partition, pedigree and recurrence are matrices
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  dimnames.recurrence <- list(c(),c("P","R"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }
  if(is.matrix(recurrence)==FALSE){
    recurrence <-  matrix(recurrence,ncol=2,dimnames=dimnames.recurrence)
  }

  new.partitions.0 <- NULL
  for(part in partition[!duplicated(partition[,"Part"]),"Part"]){

    recurrence.0 <- matrix(recurrence[recurrence[,"P"]==part,],nrow=1,ncol=2, dimnames=dimnames.recurrence)

    if(recurrence.0[,"P"]==part&recurrence.0[,"R"]==1){

      partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
      ind.i <- individual_i(pedigree, partition.subset)
      block.w.i <- matrix(partition.subset[partition.subset[,"Block"]==partition.subset[partition.subset[,"Ind"]==ind.i,"Block"][1],],ncol=5,dimnames=dimnames.partition) # block with individual i
      i.rows <- which(block.w.i[,"Ind"]==ind.i)
      block.w.i.i <- matrix(block.w.i[i.rows,],ncol=5,dimnames=dimnames.partition)
      block.w.i.others <- matrix(block.w.i[-i.rows,],ncol=5,dimnames=dimnames.partition)
      blocks.wo.i <- matrix(partition.subset[partition.subset[,"Part"]==part&partition.subset[,"Block"]!=block.w.i[,"Block"][1],],ncol=5,dimnames=dimnames.partition) # blocks without individual i

      gender.i <- pedigree[pedigree[,"iid"]==ind.i,"sex"]
      mother.i <- pedigree[pedigree[,"iid"]==ind.i,"mid"]
      father.i <- pedigree[pedigree[,"iid"]==ind.i,"pid"]


      if(gender.i==1){
        block.w.i.i.new <- matrix(block.w.i.i[1,],ncol=5,dimnames=dimnames.partition)
        block.w.i.i.new[,"Ind"] <- mother.i

        new.partitions <- rbind(block.w.i.i.new, block.w.i.others, blocks.wo.i)
        new.partitions[,"Part"] <- as.numeric(paste(part,1,sep="")) # renaming the partition
      }

      if(gender.i==2){
        no.genes <- length(block.w.i.i[,1])
        coef.1 <- 1 - 2^(1-no.genes)
        coef.2 <- 0.5^no.genes
        part.1 <- 1

        block.w.i.i.mother <- matrix(block.w.i.i[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
        block.w.i.i.mother[,"Ind"] <- mother.i
        block.w.i.i.father <- matrix(block.w.i.i[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
        block.w.i.i.father[,"Ind"] <- father.i
        block.w.i.i.mother.father <- matrix(rbind(block.w.i.i.mother,block.w.i.i.father),nrow=2,ncol=5,dimnames=dimnames.partition)

        if(no.genes!=1){
          new.partition.1 <- matrix(rbind(block.w.i.i.mother.father, block.w.i.others, blocks.wo.i),ncol=5,dimnames=dimnames.partition)
          new.partition.1[,"Coef"] <- new.partition.1[,"Coef"]*coef.1
          new.partition.1[,"Part"] <- as.numeric(paste(part,part.1,sep=""))
          part.1 <- part.1 + 1
        } else{ new.partition.1 <- NULL}

        new.partition.2 <- matrix(rbind(block.w.i.i.mother, block.w.i.others, blocks.wo.i),ncol=5,dimnames=dimnames.partition)
        new.partition.2[,"Coef"] <- new.partition.2[,"Coef"]*coef.2
        new.partition.2[,"Part"] <- as.numeric(paste(part,part.1,sep=""))
        part.1 <- part.1 + 1
        new.partition.3 <- matrix(rbind(block.w.i.i.father, block.w.i.others, blocks.wo.i),ncol=5,dimnames=dimnames.partition)
        new.partition.3[,"Coef"] <- new.partition.3[,"Coef"]*coef.2
        new.partition.3[,"Part"] <- as.numeric(paste(part,part.1,sep=""))

        new.partitions <- rbind(new.partition.1,new.partition.2,new.partition.3)
      }
      new.partitions.0 <- rbind(new.partitions.0,new.partitions)
    }
  }
  return(new.partitions.0)
}



recurrence_rule_2 <- function(pedigree, partition, recurrence){
  # Step 1. Check that partition, pedigree and recurrence are matrices
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  dimnames.recurrence <- list(c(),c("P","R"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }
  if(is.matrix(recurrence)==FALSE){
    recurrence <-  matrix(recurrence,ncol=2,dimnames=dimnames.recurrence)
  }

  new.partitions.0 <- NULL
  for(part in partition[!duplicated(partition[,"Part"]),"Part"]){

    recurrence.0 <- matrix(recurrence[recurrence[,"P"]==part,],nrow=1,ncol=2, dimnames=dimnames.recurrence)

    if(recurrence.0[,"P"]==part&recurrence.0[,"R"]==2){

      partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
      ind.i <- individual_i(pedigree, partition.subset)
      blocks.w.i <- partition.subset[partition.subset[,"Ind"]==ind.i,"Block"][!duplicated(partition.subset[partition.subset[,"Ind"]==ind.i,"Block"])]
      block.w.i.1 <- matrix(partition.subset[partition.subset[,"Block"]==blocks.w.i[1],],ncol=5,dimnames=dimnames.partition) # first block with individual i
      block.w.i.2 <- matrix(partition.subset[partition.subset[,"Block"]==blocks.w.i[2],],ncol=5,dimnames=dimnames.partition) # second block with individual i
      i.rows.1 <- which(block.w.i.1[,"Ind"]==ind.i)
      i.rows.2 <- which(block.w.i.2[,"Ind"]==ind.i)
      block.w.i.i.1 <- matrix(block.w.i.1[i.rows.1,],ncol=5,dimnames=dimnames.partition)
      block.w.i.i.2 <- matrix(block.w.i.2[i.rows.2,],ncol=5,dimnames=dimnames.partition)
      block.w.i.others.1 <- matrix(block.w.i.1[-i.rows.1,],ncol=5,dimnames=dimnames.partition)
      block.w.i.others.2 <- matrix(block.w.i.2[-i.rows.2,],ncol=5,dimnames=dimnames.partition)
      blocks.wo.i <- matrix(partition.subset[!(partition.subset[,"Block"]==blocks.w.i[1]|partition.subset[,"Block"]==blocks.w.i[2]),],ncol=5,dimnames=dimnames.partition) # blocks without individual i
      no.genes <- sum(length(block.w.i.i.1[,1]),length(block.w.i.i.2[,1]))
      coef.1 <- 0.5^no.genes

      gender.i <- pedigree[pedigree[,"iid"]==ind.i,"sex"]
      mother.i <- pedigree[pedigree[,"iid"]==ind.i,"mid"]
      father.i <- pedigree[pedigree[,"iid"]==ind.i,"pid"]

      if(gender.i==1){ stop("Recurrence rule 2 violated - i is a male")}

      block.w.i.i.mother.1 <- matrix(block.w.i.i.1[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
      block.w.i.i.mother.1[,"Ind"] <- mother.i
      block.w.i.i.father.1 <- matrix(block.w.i.i.2[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
      block.w.i.i.father.1[,"Ind"] <- father.i
      block.w.i.i.mother.2 <- matrix(block.w.i.i.2[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
      block.w.i.i.mother.2[,"Ind"] <- mother.i
      block.w.i.i.father.2 <- matrix(block.w.i.i.1[1,],nrow=1,ncol=5,dimnames=dimnames.partition)
      block.w.i.i.father.2[,"Ind"] <- father.i


      block.w.i.mother.new.1 <- rbind(block.w.i.i.mother.1,block.w.i.others.1)
      block.w.i.father.new.1 <- rbind(block.w.i.i.father.1,block.w.i.others.2)
      block.w.i.mother.new.2 <- rbind(block.w.i.i.mother.2,block.w.i.others.2)
      block.w.i.father.new.2 <- rbind(block.w.i.i.father.2,block.w.i.others.1)

      new.partition.1 <- matrix(rbind(block.w.i.mother.new.1,block.w.i.father.new.1,blocks.wo.i),ncol=5,dimnames=dimnames.partition)
      new.partition.1[,"Coef"] <- new.partition.1[,"Coef"]*coef.1
      new.partition.1[,"Part"] <- as.numeric(paste(part,1,sep=""))
      new.partition.2 <- matrix(rbind(block.w.i.father.new.2,block.w.i.mother.new.2,blocks.wo.i),ncol=5,dimnames=dimnames.partition)
      new.partition.2[,"Coef"] <- new.partition.2[,"Coef"]*coef.1
      new.partition.2[,"Part"] <- as.numeric(paste(part,2,sep=""))

      new.partitions.0 <- rbind(new.partitions.0, matrix(rbind(new.partition.1,new.partition.2),ncol=5,dimnames=dimnames.partition))
    }
  }
  return(new.partitions.0)
}




boundary_condition_1 <- function(pedigree, partition, current.phi){
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }
  B1 <- NULL
  for(part in partition[,"Part"]){
    partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
    coef <- 1
    for(ind in partition.subset[!duplicated(partition.subset[,"Ind"]),"Ind"]){
      if(pedigree[pedigree[,"iid"]==ind,"sex"]==2&sum(!duplicated(partition.subset[partition.subset[,"Ind"]==ind,"Block"]))>=3){
        coef <- 0
      }
      if(pedigree[pedigree[,"iid"]==ind,"sex"]==1&sum(!duplicated(partition.subset[partition.subset[,"Ind"]==ind,"Block"]))>=2){
        coef <- 0
      }
    }
    B1 <- rbind(B1, cbind(part, coef))
  }
  if(any(B1[,"coef"]==0)){
    new.partition <- partition[B1[,"coef"]!=0,]
    return(list(new.partition,current.phi))
  } else{ return(list(partition,current.phi))}
}


boundary_condition_2 <- function(pedigree, partition, current.phi){
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }
  founders <- pedigree[pedigree[,"pid"]==0&pedigree[,"mid"]==0,]
  B2 <- NULL
  for(part in partition[,"Part"]){
    partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
    coef <- 1
    for(block in partition.subset[,"Block"]){
      block.inds <- partition.subset[partition.subset[,"Block"]==block,"Ind"]
      block.inds <- block.inds[!duplicated(block.inds)]
      f.block <- 0
      for(ind in block.inds){
        if(any(ind==founders[,"iid"])){ f.block <- f.block + 1}
      }
      if(f.block>=2){ coef <- 0}
    }
    B2 <- rbind(B2, cbind(part, coef))
  }
  if(any(B2[,"coef"]==0)){
    new.partition <- partition[B2[,"coef"]!=0,]
    return(list(new.partition,current.phi))
  } else{ return(list(partition,current.phi))}
}


boundary_condition_3 <- function(pedigree, partition, current.phi){
  dimnames.ped <- list(c(),c("iid","pid","mid","sex"))
  dimnames.partition <- list(c(),c("Part","Block","Ind","Allele","Coef"))
  if(is.matrix(pedigree)==FALSE){
    pedigree <-  matrix(pedigree,ncol=4,dimnames=dimnames.ped)
  }
  if(is.matrix(partition)==FALSE){
    partition <-  matrix(partition,ncol=5,dimnames=dimnames.partition)
  }
  founders <- pedigree[pedigree[,"pid"]==0&pedigree[,"mid"]==0,]
  B3 <- NULL
  for(part in partition[,"Part"]){
    partition.subset <- matrix(partition[partition[,"Part"]==part,],ncol=5,dimnames=dimnames.partition)
    ind.founder.0 <- NULL
    m1 <- 0; m2.0 <- NULL
    for(ind in partition.subset[,"Ind"]){
      if(any(ind==founders[,"iid"])){
        ind.founder <- 1
        if(pedigree[pedigree[,"iid"]==ind,"sex"]==2){
          m1 <- m1 + 1
          m2.0 <- rbind(m2.0,ind)
        }
      } else{ind.founder <- 0}
      ind.founder.0 <- rbind(ind.founder.0, ind.founder)
    }
    if(all(ind.founder.0==1)){
      m2 <- length(m2.0[!duplicated(m2.0)])
      B3 <- rbind(B3, cbind(1, (1/2)^(m1 - m2)*partition.subset[1,"Coef"]))
    } else { B3 <- rbind(B3, cbind(0, 0))}
  }
  new.partition.0 <- matrix(cbind(partition, B3),ncol=7)
  new.partition.1 <- matrix(new.partition.0[new.partition.0[,6]==0,1:5],ncol=5,dimnames=dimnames.partition)
  new.partition.2 <- matrix(new.partition.0[new.partition.0[,6]==1,],ncol=7,dimnames=list(c(),c("Part","Block","Ind","Allele","Coef","b3","new.coef")))

  if(dim(new.partition.2)[1]!=0){
    new.phi <- sum(new.partition.2[!duplicated(new.partition.2[,"Part"]),7],current.phi)
  } else{ new.phi <- current.phi}
  return(list(new.partition.1,new.phi))
}



phi <- function(pedigree, partition){
  B.1 <- boundary_condition_1(pedigree, partition, 0); B.1
  B.2 <- boundary_condition_2(pedigree, B.1[[1]], B.1[[2]]); B.2
  B.3 <- boundary_condition_3(pedigree, B.2[[1]], B.2[[2]]); B.3
  B.check <- dim(B.3[[1]])[1]

  while(B.check!=0){
    R.selection <- recurrence_selection(pedigree, B.3[[1]]); R.selection
    R.1 <- recurrence_rule_1(pedigree,B.3[[1]],R.selection); R.1
    R.2 <- recurrence_rule_2(pedigree,B.3[[1]],R.selection); R.2
    R.combined <- rbind(R.1,R.2); R.combined

    B.1 <- boundary_condition_1(pedigree, R.combined, B.3[[2]]); B.1
    B.2 <- boundary_condition_2(pedigree, B.1[[1]], B.1[[2]]); B.2
    B.3 <- boundary_condition_3(pedigree, B.2[[1]], B.2[[2]]); B.3
    B.check <- dim(B.3[[1]])[1]
  }
  return(B.3[[2]])
}


psi <- function(pedigree, n, i, j){
  sex_i <- pedigree[pedigree[,"iid"] == i,"sex"]
  sex_j <- pedigree[pedigree[,"iid"] == j,"sex"]
  if(sex_i==1&sex_j==1){
    if(n==1){ return(phi(pedigree,partitions(pedigree,1,i,j))) }
    if(n==2){ return(phi(pedigree,partitions(pedigree,2,i,j))) }
  }
  if((sex_i==1&sex_j==2)|(sex_i==2&sex_j==1)){
    if(n==1){ return(phi(pedigree,partitions(pedigree,1,i,j))) }
    if(n==2){ return(phi(pedigree,partitions(pedigree,2,i,j))) }
    if(n==3){ return(2*phi(pedigree,partitions(pedigree,3,i,j))) }
    if(n==4){ return(phi(pedigree,partitions(pedigree,4,i,j))) }
  }
  if(sex_i==2&sex_j==2){
    if(n==1){ return(phi(pedigree,partitions(pedigree,1,i,j))) }
    if(n==2){ return(phi(pedigree,partitions(pedigree,2,i,j))) }
    if(n==3){ return(2*phi(pedigree,partitions(pedigree,3,i,j))) }
    if(n==4){ return(phi(pedigree,partitions(pedigree,4,i,j))) }
    if(n==5){ return(2*phi(pedigree,partitions(pedigree,5,i,j))) }
    if(n==6){ return(phi(pedigree,partitions(pedigree,6,i,j))) }
    if(n==7){ return(2*phi(pedigree,partitions(pedigree,7,i,j))) }
    if(n==8){ return(4*phi(pedigree,partitions(pedigree,8,i,j))) }
    if(n==9){ return(phi(pedigree,partitions(pedigree,9,i,j))) }
  }
}

delta <- function(pedigree, n, i, j){
  sex_i <- pedigree[pedigree[,"iid"] == i,"sex"]
  sex_j <- pedigree[pedigree[,"iid"] == j,"sex"]
  if(sex_i==1&sex_j==1){
    if(n==1){ return(psi(pedigree, 1, i, j)) }
    if(n==2){ return(psi(pedigree, 2, i, j)) } else(print("Not a valid n in calculation of delta"))
  }
  if((sex_i==1&sex_j==2)|(sex_i==2&sex_j==1)){
    if(n==1){ return(psi(pedigree, 1, i, j) - psi(pedigree, 3, i, j)/2) }
    if(n==2){ return(psi(pedigree, 2, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 4, i, j)) }
    if(n==3){ return(2*psi(pedigree, 3, i, j)) }
    if(n==4){ return(2*psi(pedigree, 4, i, j)) }  else(print("Not a valid n in calculation of delta"))
  }
  if(sex_i==2&sex_j==2){
    if(n==1){ return(psi(pedigree, 1, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 5, i, j)/2 + psi(pedigree, 7, i, j)/2 + psi(pedigree, 8, i, j)/4) }
    if(n==2){ return(psi(pedigree, 2, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 4, i, j) - psi(pedigree, 5, i, j)/2 - psi(pedigree, 6, i, j) + psi(pedigree, 7, i, j)/2 + 3*psi(pedigree, 8, i, j)/4 + psi(pedigree, 9, i, j)) }
    if(n==3){ return(2*psi(pedigree, 3, i, j) - 2*psi(pedigree, 7, i, j) - psi(pedigree, 8, i, j)) }
    if(n==4){ return(2*psi(pedigree, 4, i, j) - psi(pedigree, 8, i, j) - 2*psi(pedigree, 9, i, j)) }
    if(n==5){ return(2*psi(pedigree, 5, i, j) - 2*psi(pedigree, 7, i, j) - psi(pedigree, 8, i, j)) }
    if(n==6){ return(2*psi(pedigree, 6, i, j) - psi(pedigree, 8, i, j) - 2*psi(pedigree, 9, i, j)) }
    if(n==7){ return(4*psi(pedigree, 7, i, j)) }
    if(n==8){ return(4*psi(pedigree, 8, i, j)) }
    if(n==9){ return(4*psi(pedigree, 9, i, j)) }  else(print("Not a valid n in calculation of delta"))
  }
}

ibd_coefs_X <- function(pedigree, i, j){
  sex_i <- pedigree[pedigree[,"iid"] == i,"sex"]
  sex_j <- pedigree[pedigree[,"iid"] == j,"sex"]
  if(sex_i==1&sex_j==1){
    z0 <- delta(pedigree, 2, i, j)
    z1 <- delta(pedigree, 1, i, j)
    z2 <- 0
  }
  if((sex_i==1&sex_j==2)|(sex_i==2&sex_j==1)){
    z0 <- delta(pedigree, 2, i, j) + delta(pedigree, 4, i, j)
    z1 <- delta(pedigree, 1, i, j) + delta(pedigree, 3, i, j)
    z2 <- 0
  }
  if(sex_i==2&sex_j==2){
    z0 <- delta(pedigree, 2, i, j) + delta(pedigree, 4, i, j) + delta(pedigree, 6, i, j) + delta(pedigree, 9, i, j)
    z1 <- delta(pedigree, 3, i, j) + delta(pedigree, 5, i, j) + delta(pedigree, 8, i, j)
    z2 <- delta(pedigree, 1, i, j) + delta(pedigree, 7, i, j)
  }
  return(cbind(z0, z1, z2))
}


