# Parent-of-Origin inference for biobanks




This pipeline has been used to infer the Parent-of-Origin (PofO) of the UK Biobank individuals.

If you use (parts) of this code, please cite [Hofmeister et al., bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.11.03.467079v1)

All the summary statistics are available via our [PofO Catalog](https://www.tinyurl.com/PofOCat). If you use these summary statistics, please cite [Hofmeister et al., bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.11.03.467079v1)


For the purpose of the example, parts of this pipeline have been adapted to be run on the 1000GP chromosome 20. The resuts of this example are not usable, these are only to give an example of the input and output files. 



</br>
  &nbsp
 








#
### Step 1: If not provided, compute the relatedness using the KING software.

This step is optional as the relatedness file used for the original manuscript was available as part of the UK Biobank data.

Starting from the example variant call file (VCF, ```step1_compute_relatedness/1000GP.chr20.MAF5.vcf.gz```), this step can be run with the following command:
	
```
cd step1_compute_relatedness/
bash step0_relatedness.sh
```

KING software and documentation can be found here: https://www.kingrelatedness.com/




</br>
  &nbsp




#
### Step 2: Duos/Trios identification & surrogate parents group inference

To identify trios and duos we used pairwise kinship and IBS0 estimates up to third degree relative computed using KING and provided as part of the UK biobank study. Following Manichaikul et al. and Bycroft et al., we defined offspring-parent pairs as having a kinship coefficient between 0.1767 and 0.3535 and an IBS0 below 0.0012. We also added the condition of age difference greater than 15 years between parent-offspring pairs. We used the age and sex of the individuals to distinguish parents and offspring. For the trios, we also ensured that the two parents have different sex. We also used pairwise kinship and IBS0 estimates up to the third degree relative to identify sibling pairs (kinship between 0.1767 and 0.3535 and IBS0 above 0.0012), and second- and third-degree relatives’ pairs (kinship below 0.1767)

* For a given individual, we first identify whether parental genomes are available in the dataset using pre-defined threshold on the kinship and IBS0 estimates. 
* We then identify all its close relatives from the 2nd and 3rd degree.
* We use the relatedness in-between these 2nd-to-3rd degree relatives to cluster them into groups based on their relatedness using the igraph R package.

All these steps can be run with the following command:

```
cd step2_groups/
bash step0_run_grouping.sh
```

IGRAPH software and documentation can be found here: https://igraph.org/






</br>
  &nbsp




#
### Step 3: Groups assignments

We assigned parental status (i.e. mother or father) to groups of close relatives by examining shared IBD segments on chromosome X using XIBD, a software specifically designed to map IBD on chromosome X. This assignment was only possible for males as they inherit their only chromosome X copy from their mother: a close relative sharing IBD on chromosome X with the target is expected to be from the maternal side of the family.

As described in the original manuscript, we used a threshold at 0.1 of IBD1 coefficient to label the maternal group.

The scripts used for this analysis can be found here ```step3_xIBD/step0_run_xIBD_on_subset.sh``` but can not be run on an example file. However, examples of output files can be found here ``` step3_xIBD/Subsets/xIBD.txt ``` .


XIBD software and documentation can be found here: https://github.com/bahlolab/XIBD







</br>
  &nbsp
 



#
### Step 4: Pofo inference

* Step 4.1. IBD mapping : We designed a Hidden Markov Model (HMM), namely Target Haplotype ORigin INference (THORIN), to identify IBD sharing between the target haplotypes and a reference panel mixing haplotypes from 2 different sources: from the surrogate parents of the target (labelled as mother or father, see Step3: group assignments) and from unrelated samples. The model uses a forward-backward procedure to compute, for each allele of a target haplotype, the probability of copying the allele from (i) the surrogate mother haplotypes, (ii) the surrogate father haplotypes or (iii) unrelated haplotypes. Here, we used 100 unrelated haplotypes as decoys so that the model is not forced to systematically copy from surrogate parents. When the model copies the target haplotype from a specific surrogate parent at a given locus with high probability, we can therefore infer the PofO at this locus from the parental group the surrogate parent belongs to. When the model copies from unrelated haplotypes, no inference can be made at the locus.


* Step 4.2. scaffold construction : We built a haplotype scaffold comprising all alleles for which we know the PofO from IBD sharing with surrogate parents. As described in the original manuscript, we only included in the scaffold IBD tracks longer than 3cM.


All these steps can be run with the following command:

```
cd step4_pofo_inference
bash step0_PofO_inference.sh
bash step1_scaffold.sh
```

THORIN software and documentation can be found in ```step4_pofo_inference/THORIN```






</br>
  &nbsp




#
### Step 5: Extrapolation by phasing

We use the scaffold previously built to rephase the data using SHAPEIT4. The goal of this second round of phasing is : (i) to ensure that the pool of alleles coming from the same parent land onto the same haplotype, (ii) to propagate the PofO assignment from IBD tracks to all alleles along the chromosomes and (ii) to correct long range switch errors. Point (ii) is made possible as all alleles with PofO unknown (i.e. not in IBD tracks) are phased relatively to the haplotype scaffold so that we can extrapolate their PofO from the scaffold they co-localize with (paternal/maternal).

In practice, we first rephase the data using the scaffold to produce a temporary tree storing all the possible phasing possibilities (step 5.1). Then, we use this tree to randomly sample 1,000 different haplotypes (step 5.2). We use the step 5.2 to identify, over the 1,000 phasing, haplotype segment being phased >=70% on the same haplotype, that we defined as "high confidence segments". Finally, we use this tree to output the most likely haplotype, as in a usual run a phasing (step 5.4).
	

All these steps can be run with the following command:
```
cd step5_pofo_confidence_score
bash step1_phasing_confidence.sh
```

SHAPEIT4 software and documentation can be found here: https://odelaneau.github.io/shapeit4/








</br>
  &nbsp
 



#
### Step 6: Extrapolation by imputation

We inferred PofO for untyped alleles, i.e. not included on the SNP array. To do so, we imputed the data using IMPUTE5 v1.1.4 with the Haplotype Reference Consortium as a reference panel. As our data is phased with each haplotype assigned to a specific parent, we used the parameter --out-ap-field to run a haploid imputation of the data and separately imputed the paternal haplotype and the maternal haplotype. As a result of haploid imputation, the PofO of imputed alleles can be probabilistically deduced from the imputation dosages: an allele imputed with a dosage of 0.85 on the paternal haplotype has 85% probability of being inherited from the father.

IMPUTE5 software and documentation can be found here : https://jmarchini.org/software/#impute-5 

</br>
  &nbsp
  </br>
  &nbsp





Software versions used in the original manuscript:

* BCFtools v1.8
* SHAPEIT v4.2.1
* BOLT-LMM v2.3.4
* PINK v1.90b5
* R v3.5.1
* R IGRAPH v1.2.2
* IMPUTE5 v1.1.4





