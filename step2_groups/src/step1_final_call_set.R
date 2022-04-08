DATA='UKB'
DATA='1000GP'

## Step0. Association set
ped<-as.data.frame(data.table::fread('Grouping/PED.txt', hea=F))
full_ped<-ped
trios<-as.data.frame(data.table::fread('Grouping/Trios.txt',hea=F))
duos<-as.data.frame(data.table::fread('Grouping/Duos.txt',hea=F))
duos$V3<-''

'%ni%'=Negate('%in%')
sub<-ped[ped$V1 %ni% trios$V1,]
sub<-sub[sub$V1 %ni% duos$V1,]
ped<-rbind(sub, trios, duos)

write.table(ped,'Grouping/Full_callset.txt', quote=F, col.names=F, row.names=F, sep='\t')



#related:
g1<-unlist(lapply(ped$V2, FUN=function(x){ unlist(strsplit(x,'='))[2]    }))
g1.2<-unlist(lapply(g1, FUN=function(x){ unlist(strsplit(x,';'))    }))
g2<-unlist(lapply(ped$V3, FUN=function(x){ unlist(strsplit(x,'='))[2]    }))
g2.2<-unlist(lapply(g2, FUN=function(x){ unlist(strsplit(x,';'))    }))
samples<-unique(c(ped$V1, g1.2, g2.2))
write.table(data.frame(V1=samples), 'Grouping/Related_samples.txt', quote=F, col.names=F, row.names=F)




## unrelated:
relative_file='../step1_compute_relatedness/King/chr1.kin0'

r<-read.table(relative_file, hea=T)
r$id1<-paste0(r$ID1,'_',r$ID1)
r$id2<-paste0(r$ID2,'_',r$ID2)
'%ni%'=Negate('%in%')

if (DATA=='1000GP'){
r<-r[r$InfType!='UN',]
related_samples<-unique(c(r$id1, r$id2))
sample_file='../step1_compute_relatedness/Plink/chr1.fam'
samples<-read.table(sample_file, hea=F); samples$V1<-paste0(samples$V1,'_', samples$V1)
unrelated<-samples$V1[samples$V1 %ni% related_samples]
unrelated<-unlist(lapply(unrelated, FUN=function(x){ unlist(strsplit(x,'_'))[1]  }))
}

if (DATA=='UKB'){
related_samples<-unique(c(r$id1, r$id2))
ukb_british_samples='' # file with british and irish sample ids
b<-read.table(ukb_british_samples, hea=F)
unrelated<-b[b$V1 %ni% related_samples,]
}


write.table(data.frame(V1=sample(unrelated,100)), 'Grouping/100unrelated.txt', quote=F, col.names=F, row.names=F)














