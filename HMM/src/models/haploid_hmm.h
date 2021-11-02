/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/
#ifndef _HAPLOID_HMM_H
#define _HAPLOID_HMM_H

#include <utils/otools.h>
#include <objects/hmm_parameters.h>
#include <containers/haplotype_set.h>

class haploid_hmm {
private:
	//EXTERNAL DATA
	hmm_parameters & M;
	haplotype_set & H;

	int ngroups;
	int idx_tar_vcf;	//Index of target haplotype in VCF
	int idx_tar_grp;	//Index of target haplotype in Groups
	vector < pair < int, int > > idxH;

	//DYNAMIC ARRAYS
	vector < float > AlphaSum;
	vector < vector < float > > Alpha;
	vector < float > Beta;

public:
	//CONSTRUCTOR/DESTRUCTOR
	haploid_hmm(int , int, haplotype_set &, hmm_parameters &);
	~haploid_hmm();

	double forward();
	double backward();
};


#endif
