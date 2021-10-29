Parent-of-origin inference from close relatives




Usage:

		$ bin/pooer -I examples/1000G.chr20.related.vcf.gz -H examples/1000G.chr20.unrelated.vcf.gz -M maps/chr20.b37.gmap.gz -G examples/group_file.txt -R 20 -O examples/output_test.prob.gz




-I : vcf file with all individuals included in the group file

-H : vcf file with individuals unrelated to all the individuals included in -I

-M : genetic map

-G : group file (col1= target sample; col2= first group of relatives; col3= second group of relatives, if any)

-R : region

-O : output file




