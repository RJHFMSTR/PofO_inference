import gzip
import os

DATA='UKB'
DATA='1000GP'
male_targets='' # Not used for the example on the 1000GP data. Used for the inference on the UKB samples.
group_file='../step2_groups/Grouping/PED.txt' # group file with trios/duos having group of relatives to allow the benchmark of the IBD X estimate.
ukb_chrX=''


## 1. get males individuals

if DATA=='UKB':
	males=dict()
	with gzip.open(male_targets,'rb') as f:
		for line in f:
			tmp=line.split()
			males[tmp[0]]=0



## 2. get male target and their relatives

	d2=dict()

	d=dict()
	file=open(group_file)
	for line in file:
		tmp=line.split()
		if tmp[0] in males:
			
		        d[tmp[0]]=[[],[]]
        
	                t1=tmp[1].split(';')
	                t11=t1[0].split('=')[1]
	                d[tmp[0]][0].append(t11)
		
	                if len(t1)>1:
	                        
	                        for i in t1[1:]:
	                                d[tmp[0]][0].append(i)
	                                
       		        if len(tmp)>2:
       	                	t2=tmp[2].split(';')
       	                	t22=t2[0].split('=')[1]
       	                	d[tmp[0]][1].append(t22)

                        	if len(t2)>1:
                        	        for i in t2[1:]:
                        	                d[tmp[0]][1].append(i)
	file.close()	
	

if DATA=='1000GP':


	count=0
## 2. Let's assume they are all males for the purpose of the example

	d2=dict()

	d=dict()
	file=open(group_file)
	for line in file:
		count+=1

		if count>0: # run XIBD only for a subset of the target as an example
			tmp=line.split()
	
			d[tmp[0]]=[[],[]]
	
			t1=tmp[1].split(';')
			t11=t1[0].split('=')[1]
			d[tmp[0]][0].append(t11)
	
			if len(t1)>1:
	
				for i in t1[1:]:
					d[tmp[0]][0].append(i)
	
			if len(tmp)>2:
				t2=tmp[2].split(';')
				t22=t2[0].split('=')[1]
				d[tmp[0]][1].append(t22)
	
				if len(t2)>1:
					for i in t2[1:]:
						d[tmp[0]][1].append(i)
	file.close()





os.system('mkdir -p Subsets')


outfile=open('Subsets/Pairs.txt','w')
outfile2=open('Subsets/Samples.txt','w')



d2=dict()
for k,v in d.items():
	g1=v[0]
	g2=v[1]
	
	#	outfile2.write(k.split('_')[0]+'\t'+k.split('_')[0]+'\n')
	d2[k.split('_')[0]]=0
	for i in g1:
		outfile.write(k+'\t'+i+'\t'+'G1'+'\n')
		#outfile2.write(i.split('_')[0]+'\t'+i.split('_')[0]+'\n')
		d2[i.split('_')[0]]=0
	if len(g2)!=0:
		for i in g2:
			outfile.write(k+'\t'+i+'\t'+'G2'+'\n')
			#outfile2.write(i.split('_')[0]+'\t'+i.split('_')[0]+'\n')
			d2[i.split('_')[0]]=0
outfile.close()



for k,v in d2.items():
	outfile2.write(k+'\t'+k+'\n')

outfile2.close()

##  4. convet the files

os.system('plink --vcf '+ukb_chrX+' --keep Subsets/Samples.txt --recode 12 --out Subsets/Plink.chrX')



## run XIBD software

#os.system('Rscript src/run_Xibd.R')











