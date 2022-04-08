import gzip
import argparse
import os


parser=argparse.ArgumentParser()
parser.add_argument("-chr","--chromosome", help="choose chr", type=int)
args=parser.parse_args()

CHR=str(args.chromosome)



group_file='step2_groups/Grouping/Full_callset.txt'
scaffolded_sampled='../step4_pofo_inference/Scaffold/chr'+CHR+'.scaffolded_samples.txt'




# 1. get all scaffolded samples
scaff=dict()
file=open(scaffolded_sampled)
for line in file:
	tmp=line.split()
	scaff[tmp[0]]=0
file.close()
print(len(scaff))

#2. get the target of interest (i.e, trios only, trios + duos, trios + duos + xtargets)


c=0
Targets=[]

file=open(group_file)
for line in file:
	c+=1
	tmp=line.split()
	if start<=c<=end:
		if tmp[0] in scaff:
			Targets.append(tmp[0])
file.close()
print(len(Targets))

os.system('mkdir -p Samples_files/')
outfile=open('Samples_files/chr'+CHR+'.txt','w')
for t in Targets:
	outfile.write(t+'\n')
outfile.close()








print('Extracting segments .....')
d={}

sampling_file='Sampling/Bin.sampling.chr'+CHR+'.vcf.gz'
with gzip.open(sampling_file) as f:
	for line in f:
		if line.find('##')==0:
			pass
		elif line.find('#CHROM')==0:
			ss=line.split()
			for t in Targets:
				if t in ss:
					d[t]=[[],[]]
				else:
					d[t]='NA'
		else:
			tmp=line.split()
			pos=tmp[1]
			
			for target in Targets:
								
				if target in ss:
					prob=tmp[ss.index(target)].split(',')
					if (prob[0]!='1' and prob[3]!='1'): # het site
						d[target][0].append(int(pos))
						d[target][1].append(max([float(prob[1]), float(prob[2])]))
						






print('Building consecutive haplotype segment...')			
os.system('mkdir -p Segments/chr'+CHR)
outfile=open('/Segments/chr'+CHR+'/chr'+CHR+'.txt','w')
for target,v in d.items():
	position=v[0]; prob=v[1]
	outfile1=open('junk.txt','w')
	for k in range(len(position)):
		outfile1.write(str(position[k])+'\t'+str(prob[k])+'\n')
	outfile1.close()
	S=[]; seg=[]
	for i in range(len(prob)):
		p=prob[i]
		if p<0.7:
#			print(p);print(prob.index(p))
			pos=position[i]
			pos_minus=pos
			pos_plus=pos
			if i!=0:
				pos_minus=position[i-1]
			if i!=(len(prob)-1):
				pos_plus=position[i+1]
			seg.append(pos); seg.append(pos_minus); seg.append(pos_plus)

			if pos==position[len(position)-1]:
				if len(seg)!=0:
					start_seg=min(seg); end_seg=max(seg)
					S.append(str(start_seg)+'-'+str(end_seg)+'/'+str(len(set(seg))))
					seg=[]	

		else:
			if len(seg)!=0:
				start_seg=min(seg); end_seg=max(seg)
				S.append(str(start_seg)+'-'+str(end_seg)+'/'+str(len(set(seg))))
				seg=[]

	w=[target]+[','.join(S)]
	outfile.write('\t'.join(w)+'\n')

outfile.close()

os.system('bgzip -f Segments/chr'+CHR+'/chr'+CHR+'.txt')

