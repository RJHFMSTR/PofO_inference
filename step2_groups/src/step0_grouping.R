`%ni%` = Negate(`%in%`)
library(igraph)


args = commandArgs(trailingOnly=TRUE)

DATA='UKB'
DATA='1000GP'


### 0. Initialization, reading input files, building initial subsets

# 0.0 Input files
relatives_file='../step1_compute_relatedness/King/chr20.kin0' # output of KING relatedness software

if (DATA=='UKB'){
samples_to_consider='' # list all samples that we are interested in for the PofO inference. In UKB, we focused only on British and Irish samples for example
age_file='' # Not used for the example on the 100GP data. Used for the PofO inference in the UK biobank. File format: col1=sample_id; col2=age
sex_file='' # Not used for the example on the 100GP data. Used for the PofO inference in the UK biobank. List only male individuals: col1=male sample_id
}


# 0.1 Relatedness file: read relatedness file and fix sample id if not consistent between genotype file (VCF/BCF) and relatedness file

rel<-as.data.frame(data.table::fread(relatives_file, hea=T))

if (DATA=='1000GP'){rel<-rel[rel$InfType!='UN',];rel<-rel[,c(1,3,6,7,10)]}

names(rel)<-c('ID1','ID2','HetHet','IBS0','Kinship')

if (DATA=='UKB'){rel$FID1<-paste0(rel$ID1,'_',rel$ID1) ; rel$FID2<-paste0(rel$ID2,'_',rel$ID2) }
if (DATA=='1000GP'){rel$FID1<-rel$ID1 ; rel$FID2<-rel$ID2}


# 0.2 Samples of interest: read sample file and keep only those one in the relatedness file. Used only for the UKB dataset.

if (DATA=='UKB'){
target_samples<-as.data.frame(data.table::fread(samples_to_consider, hea=F))
rel<-rel[rel$FID1 %in% target_samples$V1 & rel$FID2 %in% target_samples$V1,]
}

samples<-unique(c(unique(rel$FID1), unique(rel$FID2)))
rel<-rel[rel$Kinship<0.3535, ] # remove MZ twins

# 0.3 Age:
if (DATA=='UKB'){
age<-as.data.frame(data.table::fread(age_file, hea=F))
names(age)<-c('id','age'); age$age<-as.numeric(age$age)}

# 0.4 Sex:
if (DATA=='UKB'){
sex<-as.data.frame(data.table::fread(sex_file, hea=F))}



# 0.5 Initialize output matrixes
Duos<-matrix(ncol=3)
Trios<-matrix(ncol=3)
Sibs<-matrix(ncol=3)
OUT=matrix(ncol=3)



### 1. Grouping of relatives
print(length(samples))
c=0
for (t in samples){
	c=c+1
	print(c)
	sex_target='NA'
	if (DATA=='UKB'){sex_target='Female'; if (t %in% sex$V1) {sex_target='Male'}}
	
	# 1.1 Grouping of first-degree relatives (parents, siblings)

	# 1.1.1 Identify parents (PO)
	PO<-unique(c(as.character(rel$FID2[rel$FID1==t & rel$Kinship>=0.1767 & rel$IBS0<0.0012]),as.character(rel$FID1[rel$FID2==t & rel$Kinship>=0.1767 & rel$IBS0<0.0012])))

	# 1.1.1.1 Duos	
	if (length(PO)==1){

		if (DATA=='UKB'){ # verify difference in age
		age1<-age$age[age$id==PO[1]]
		if (age$age[age$id==t]+15<age1){
			parent1=PO[1]
			sex1='Female'
			if (PO[1] %in% sex$V1) {sex1='Male'}
			duos<-matrix(c(t, paste0('G1=',PO[1] ), sex1), ncol=3)
			Duos<-rbind(Duos, duos) }
			}

		if (DATA=='1000GP'){
			parent1=PO[1]
                        sex1='NA'
                        duos<-matrix(c(t, paste0('G1=',PO[1] ), sex1), ncol=3)
                        Duos<-rbind(Duos, duos)	
			}
		}


  	# 1.1.1.2 Trios
	if (length(PO)==2){
		
		if (DATA=='UKB'){
		age1<-age$age[age$id==PO[1]]
		age2<-age$age[age$id==PO[2]]
		if ((age$age[age$id==t]+15<age1) & (age$age[age$id==t]+15<age2)){
			parent1=PO[1]
			sex1='Female'
			sex2='Female'
			if (PO[1] %in% sex$V1) {sex1='Male'}
			if (PO[2] %in% sex$V1) {sex2='Male'}
			if (sex1!=sex2){
				Tmp<-data.frame(id=c(PO[1], PO[2]), sex=c(sex1, sex2))
				trios<-matrix(c(t, paste0('G1=', Tmp$id[Tmp$sex=='Male']), paste0('G2=', Tmp$id[Tmp$sex=='Female'])), ncol=3)
				Trios<-rbind(Trios, trios)  }}
			}
		
		if(DATA=='1000GP'){
                        sex1='NA'
                        sex2='NA'
			trios<-matrix(c(t, paste0('G1=', PO[1]), paste0('G2=', PO[2])), ncol=3)
			Trios<-rbind(Trios, trios)
                	}	
  		}	




	# 1.1.2 Identify Siblings (FS)
	FS<-unique(c(as.character(rel$FID2[rel$FID1==t & rel$Kinship>=0.1767 & rel$IBS0>=0.0012]),as.character(rel$FID1[rel$FID2==t & rel$Kinship>=0.175 & rel$IBS0>=0.0012])))
	sibs<-matrix(c(t, paste(FS, collapse=';'), sex_target), ncol=3)
        Sibs<-rbind(Sibs, sibs)	


	# 1.2 Grouping of second- and third-degree relatives

	g1<-""; g2<-""
	relatives<-unique(c(as.character(rel$FID2[rel$FID1==t & rel$Kinship<0.1767]),as.character(rel$FID1[rel$FID2==t & rel$Kinship<0.1767])))
	
	if (length(relatives)!=0){

		d_pairs<-data.frame()
		for (relat in relatives){
			relatives2<-unique(c(as.character(rel$FID2[rel$FID1==relat]),as.character(rel$FID1[rel$FID2==relat])))
			relatives2<-relatives2[relatives2 %in% relatives]; relatives2<-relatives2[relatives2!=t]
			d_pairs<-rbind(d_pairs, data.frame(t1=rep(relat, length(relatives2)), t2=relatives2))
			if (length(relatives2)==0){d_pairs<-rbind(d_pairs, data.frame(t1=relat, t2=relat)) }}   

		m<-as.matrix(d_pairs); g<-simplify(graph_(m, from_edgelist(), directed = FALSE));
		dg <- decompose.graph(g)
    
		p='T'
		if (length(dg)==0){p='F'}
		if (length(dg)==1) { g1=vertex_attr(dg[[1]])$name; g2=""}
		if (length(dg)==2){ g1=vertex_attr(dg[[1]])$name; g2=vertex_attr(dg[[2]])$name}
		if (length(dg)>2){p='F'}
    
		if (g1[1]=="" & g2[1]==""){ p='F'}
		if (g1[1]!="" & g2[1]==""){ G1<-paste0("G1=",paste0(g1, collapse=';')); G2=""; out=matrix(c(t,G1, G2), ncol=3); OUT=rbind(OUT, out)}
		if (g1[1]!="" & g2[1]!="") { G1<-paste0("G1=",paste0(g1, collapse=';'));  G2<-paste0("G2=",paste0(g2, collapse=';')); out=matrix(c(t,G1, G2), ncol=3); OUT=rbind(OUT, out)}
    
		d_pairs2<-d_pairs
		for (relat in relatives){d_pairs2<-rbind(d_pairs2, data.frame(t1=t, t2=relat))}

		m2<-as.matrix(d_pairs2); g2<-simplify(graph_(m2, from_edgelist(), directed = FALSE))
		V(g2)$color[names(V(g2))==t]<-'green'

#		png(paste0(t,'_grouping.png'))
#		plot(g2)
#		dev.off()
	}
   
}






### 2. Output files:

OUT<-matrix(OUT[-1, ], ncol=3)
Sibs<-matrix(Sibs[-1, ], ncol=3)
Duos<-matrix(Duos[-1, ], ncol=3)
Trios<-matrix(Trios[-1, ], ncol=3)


system('mkdir -p Grouping')
write.table(OUT, paste0("Grouping/PED.txt"), quote=F, row.names=F, col.names=F, sep='\t')
write.table(Trios, paste0("Grouping/Trios.txt"), quote=F, row.names=F, col.names=F, sep='\t')
write.table(Duos, paste0("Grouping/Duos.txt"), quote=F, row.names=F, col.names=F, sep='\t')
write.table(Sibs, paste0("Grouping/Sibs.txt"), quote=F, row.names=F, col.names=F, sep='\t')





