`%ni%` = Negate(`%in%`)
# load XIBD library
library(identity, lib.loc='./')
library(XIBD, lib.loc='./')





ped<-as.data.frame(data.table::fread("Subsets/Plink.chrX.ped", hea=F))
map<-as.data.frame(data.table::fread('Subsets/Plink.chrX.map', hea=F))
pedmap<-list(ped,map)
my_genotypes <- getGenotypes(ped.map = pedmap,
                             model = 1)


# estimate parameters
my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, 
                                  number.cores = 1)


d_ibd_parent<-my_parameters[[2]]
d_ibd_parent$id1<-paste0(d_ibd_parent$fid1,'_', d_ibd_parent$fid1)
d_ibd_parent$id2<-paste0(d_ibd_parent$fid2,'_', d_ibd_parent$fid2)


write.table(d_ibd_parent, "Subsets/xIBD.txt", quote=F, col.names=T, row.names=F)

