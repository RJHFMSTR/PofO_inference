import argparse
import gzip
import os
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-chr','--chromosome', type=int)
parser.add_argument('-cm','--centimorgan', type=float)

args = parser.parse_args()


CHR=str(args.chromosome)
MIN=args.centimorgan

print('Initialization: processing chromosome '+CHR)


DATA='UKB'
DATA='1000GP'


Group_file='../step2_groups/Grouping/Full_callset.txt'
Prob_file='Out/chr'+CHR+'/chr'+CHR+'.prob.gz'
Scaffold_in='Input/chr'+CHR+'/related.chr'+CHR+'.vcf.gz'
Scaffold_out='Scaffold/chr'+CHR+'/chr'+CHR+'.scaffold.vcf'
scaffolded_samples='Scaffold/chr'+CHR+'/chr'+CHR+'.scaffolded_samples.txt'
os.system('mkdir -p Scaffold/chr'+CHR)



# 1. get nbr of groups for each target
d=dict()
dd=dict()
file=open(Group_file)
for line in file:
	if line.find('Target')==0:
		pass
	else:
		tmp=line.split()
		d[tmp[0]]=len(tmp)
		dd[tmp[0]]=tmp
file.close()



# 2. read and store PofO probabilities
# 2.0 initialize output dict
targets=d.keys()


if DATA=='1000GP':
	targets=targets[:10]

dt1=dict()
for i in targets:
	dt1[i]=dict()

target2=[]
pos=[]
d_pos=dict()


# 2.1 compute joint h0/h1 probabilities
print('Computing probabilities....')
with gzip.open(Prob_file,'rb') as f:
	for line in f:
		tmp=line.split()
		if line.find('#')==0:
			ss=tmp
		else:

			p=int(tmp[1])
			pos.append(p)
			d_pos[p]=float(tmp[3])
			w=tmp[:4]
			for t in targets:
				t2=t+'_HOLE_0'
				if t2 in ss:
					target2.append(t)
					nbr_group=d[t]
					if nbr_group==2:
						h1g1=float(tmp[ss.index(t+'_G1_0')])
						h1U=float(tmp[ss.index(t+'_HOLE_0')])
	                                        h2g1=float(tmp[ss.index(t+'_G1_1')])
	                                        h2U=float(tmp[ss.index(t+'_HOLE_1')])
	                                        
						h1g2=0
						h2g2=0
	                                if nbr_group==3:
	                                        h1g1=float(tmp[ss.index(t+'_G1_0')])
       		                                h1U=float(tmp[ss.index(t+'_HOLE_0')])
       	                                	h2g1=float(tmp[ss.index(t+'_G1_1')])
                                        	h2U=float(tmp[ss.index(t+'_HOLE_1')])
                                        	h1g2=float(tmp[ss.index(t+'_G2_0')])
                                        	h2g2=float(tmp[ss.index(t+'_G2_1')])
					PA=float(h1g1)*float(h2g2)+float(h1g1)*float(h2U)+float(h1U)*float(h2g2)
					PB=float(h1g2)*float(h2g1)+float(h1g2)*float(h2U)+float(h1U)*float(h2g1)
					PC=float(h1U)*float(h2U)
					PD=float(h1g1)*float(h2g1)+float(h1g2)*float(h2g2)
                                	PA2=PA/(PA+PB+PC+PD)
                                	PB2=PB/(PA+PB+PC+PD)
                                	PC2=PC/(PA+PB+PC+PD)
                                	PD2=PD/(PA+PB+PC+PD)
					P=[PA2,PB2,PC2,PD2]
					dt1[t][p]=[max(P), P.index(max(P))] # store the joint prob for each pos

# 3. identify consecutive segment of the same PofO over 3cM or more
print('Building scaffold ....')
target=target2


pos.sort()
d_scaffold=dict()

for t in targets:
	d=dt1[t];d_scaffold[t]=dict();jp=[];a=[];
	for i in pos:
		jp.append(d[i][0]);a.append(d[i][1])
	max_iter=len(a)-1;init_index=0;init_prob=a[0];jp_count=0;accepted=[0,1];max_small_prob=5;scaffold_nbr=1
	for i in range(1, len(a)):
		sub_a=a[init_index:i+1]
		if jp[i]<0.9:
			jp_count+=1	
		if (sub_a.count(init_prob)>=len(sub_a)-1) and jp_count>=max_small_prob:
			pass
		if (sub_a.count(init_prob)<len(sub_a)-1) or jp_count>=max_small_prob:
			good_sub_a=a[init_index:i]
			end_index=i-1
	                if init_prob in accepted:
	                        if (d_pos[pos[end_index]]-d_pos[pos[init_index]])>=MIN:
					index_pos_to_include=[index for index, value in enumerate(good_sub_a) if value == init_prob]
					for j in range(len(pos[init_index: end_index+1])):
						if j in index_pos_to_include:
							d_scaffold[t][pos[init_index: end_index+1][j]]=init_prob
					scaffold_nbr+=1
			init_prob=a[i];	init_index=i;jp_count=0
	        elif i==max_iter:
			end_index=i     
			if init_prob in accepted:
				if (d_pos[pos[end_index]]-d_pos[pos[init_index]])>=MIN:
					index_pos_to_include=[index for index, value in enumerate(sub_a) if value == init_prob]
			                for j in range(len(pos[init_index: end_index+1])):
			                        if j in index_pos_to_include:
				                        w=[str(pos[init_index: end_index+1][j]), str(init_prob), str(scaffold_nbr)]
							d_scaffold[t][pos[init_index: end_index+1][j]]=init_prob



# 4. write scaffold file
print('Writing scaffold...')	
D=dict()
outfile=open(Scaffold_out,'w')
with gzip.open(Scaffold_in,'rb') as f:
	for line in f:
		tmp=line.split()
		if line.find('#')==0:
			if line.find('#CHROM')==0:
				ss=tmp
				w=tmp[:9]
				for t in targets:
					w.append(t)
				outfile.write('\t'.join(w)+'\n')
			else:
				outfile.write(line)
		else:
			pos=int(tmp[1]);w=tmp[:9]
			for t in targets:
				index_target=ss.index(t);geno_target=tmp[index_target];	dt=d_scaffold[t]
				if pos in dt:
					if int(geno_target.split('|')[0])+int(geno_target.split('|')[1])==1:
						g=geno_target.split('|')[0]+'/'+geno_target.split('|')[1]
						if dt[pos]==0:
                                                        g=geno_target.split('|')[1]+'|'+geno_target.split('|')[0]
                                                        D[t]=0
						elif dt[pos]==1:
                                                        g=geno_target.split('|')[0]+'|'+geno_target.split('|')[1]
                                                        D[t]=0
					else:
						g=geno_target.split('|')[0]+'/'+geno_target.split('|')[1]
				else:
					g=geno_target.split('|')[0]+'/'+geno_target.split('|')[1]

				w.append(g)
			outfile.write('\t'.join(w)+'\n')
outfile.close()
os.system('bgzip -f '+Scaffold_out)
os.system('tabix -p vcf -f '+Scaffold_out+'.gz')



outfile=open(scaffolded_samples,'w')
for k,v in D.items():
	outfile.write(k+'\n')
outfile.close()

os.system('bgzip -f '+scaffolded_samples)
os.system('rm '+Scaffold_in)
