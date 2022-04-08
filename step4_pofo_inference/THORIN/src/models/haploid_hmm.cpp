#include <models/haploid_hmm.h>

haploid_hmm::haploid_hmm(int s, int p, haplotype_set & _H, hmm_parameters & _M) : H(_H), M(_M) {
	//Indexing of the target haplotype
	idx_tar_grp = 2 * s + p;
	idx_tar_vcf = 2 * H.targetIDXs[s] + p;

	//Push in idxH conditioning haplotype indexes coming from the groups
	idxH.clear();
	ngroups = H.sourceIDXs[s].size();
	for (int g = 0 ; g < ngroups ; g ++) {
		for (int h = 0 ; h < H.sourceIDXs[s][g].size() ; h ++) {

			int h0 = H.sourceIDXs[s][g][h]*2+0;
			int h1 = H.sourceIDXs[s][g][h]*2+1;
			assert(h0 >= 0 && h0 < H.n_hap);
			assert(h1 >= 0 && h1 < H.n_hap);
			assert(h0 != idx_tar_vcf);
			assert(h1 != idx_tar_vcf);

			idxH.push_back(pair < int, int > (H.sourceIDXs[s][g][h]*2+0, g));
			idxH.push_back(pair < int, int > (H.sourceIDXs[s][g][h]*2+1, g));
		}
	}

	//Push in idxH conditioning haplotype indexes coming from the hole/reference
	ngroups += (H.IDs_map.size()*2 != H.n_hap);
	for (int h = H.IDs_map.size() * 2 ; h < H.n_hap && ngroups;h ++)
		idxH.push_back(pair < int , int > (h, ngroups-1));

	//Sort idxH given first field (i.e. haplotype index in bitmatrix
	sort(idxH.begin(), idxH.end());

	//Allocate vectors for forward/backward computations
	Alpha = vector < vector < float > > (H.n_site, vector < float > (idxH.size(), 0.0));
	AlphaSum = vector < float > (H.n_site, 0.0);
	Beta = vector < float > (idxH.size(), 1.0);
}

haploid_hmm::~haploid_hmm() {
	idx_tar_grp = -1;
	idx_tar_vcf = -1;
	idxH.clear();
	Alpha.clear();
	AlphaSum.clear();
	Beta.clear();
}

double haploid_hmm::forward() {
	double loglik = 0.0, fact1, fact2;
	for (int l = 0 ; l < H.n_site ; l ++) {
		AlphaSum[l] = 0.0;
		bool a = H.H_opt_var.get(l, idx_tar_vcf);
		if (l == 0) {
			fact1 = 1.0 / idxH.size();
			for (int k = 0 ; k < idxH.size() ; k ++) {
				float emit = (a==H.H_opt_var.get(l,idxH[k].first))?0.9999:0.0001;
				Alpha[l][k] = emit * fact1;
				AlphaSum[l] += Alpha[l][k];
			}
		} else {
			fact1 = M.t[l-1] / idxH.size();
			fact2 = M.nt[l-1] / AlphaSum[l-1];
			for (int k = 0 ; k < idxH.size() ; k ++) {
				float emit = (a==H.H_opt_var.get(l, idxH[k].first))?0.9999:0.0001;
				Alpha[l][k] = (Alpha[l-1][k] * fact2 + fact1) * emit;
				AlphaSum[l] += Alpha[l][k];
			}
		}
		loglik += log (AlphaSum[l]);
	}
	return loglik;
}

double haploid_hmm::backward() {
	double loglik = 0.0, fact1, fact2, betaSumPrev, betaSumCurr;
	for (int l = H.n_site-1 ; l >= 0 ; l --) {
		betaSumCurr = 0.0;
		bool a = H.H_opt_var.get(l, idx_tar_vcf);
		if (l == H.n_site - 1) {
			for (int k = 0 ; k < idxH.size() ; k ++) {
				float emit = (a==H.H_opt_var.get(l,idxH[k].first))?0.9999:0.0001;
				betaSumCurr += emit;
			}
		} else {
			fact1 = M.nt[l] / betaSumPrev;
			for (int k = 0 ; k < idxH.size() ; k ++) {
				float curr_emit = (a==H.H_opt_var.get(l,idxH[k].first))?0.9999:0.0001;
				float next_emit = (H.H_opt_var.get(l+1,idx_tar_vcf)==H.H_opt_var.get(l+1,idxH[k].first))?0.9999:0.0001;
				Beta[k] = Beta[k] * fact1 * next_emit + M.t[l];
				betaSumCurr += curr_emit * Beta[k];
			}
		}

		// Expectation pass
		float sum = 0.0;
		for (int k = 0 ; k < idxH.size() ; k ++) {
			double val = Alpha[l][k] * Beta[k];
			H.copyingProbabilities[idx_tar_grp][ngroups * l + idxH[k].second] += val;
			sum += val;
		}
		for (int g = 0 ; g < ngroups ; g ++) H.copyingProbabilities[idx_tar_grp][ngroups * l + g] /= sum;

		betaSumPrev = betaSumCurr / idxH.size();
		loglik += log (betaSumPrev);
	}
	return loglik;
}
