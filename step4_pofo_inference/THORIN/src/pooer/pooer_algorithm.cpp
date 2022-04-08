////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <pooer/pooer_header.h>

#include <models/haploid_hmm.h>

void * pooe_callback(void * ptr) {
	pooer * S = static_cast< pooer * >( ptr );
	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->H.targetIDXs.size()) S->pooe(id_job);
		else pthread_exit(NULL);
	}
}

void pooer::pooe (int id_job) {
	vrb.title("Processing [" + H.IDs[H.targetIDXs[id_job]] + "]");
	vrb.bullet("Create HMM for hap0");
	haploid_hmm HM0(id_job, 0, H, M);
	vrb.bullet("Forward HMM for hap0 ...");
	double loglikelihood = HM0.forward();
	vrb.bullet("\tLLfor="+stb.str(loglikelihood, 3));
	vrb.bullet("Backward HMM for hap0");
	loglikelihood = HM0.backward();
	vrb.bullet("\tLLbac="+stb.str(loglikelihood, 3));
	vrb.bullet("Create HMM for hap1");
	haploid_hmm HM1(id_job, 1, H, M);
	vrb.bullet("Forward HMM for hap1");
	loglikelihood = HM1.forward();
	vrb.bullet("\tLLfor="+stb.str(loglikelihood, 3));
	vrb.bullet("Backward HMM for hap1");
	loglikelihood = HM1.backward();
	vrb.bullet("\tLLbac="+stb.str(loglikelihood, 3));
}

void pooer::pooe() {
	tac.clock();
	int n_thread = options["thread"].as < int > ();
	i_workers = 0; i_jobs = 0;
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, pooe_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < H.targetIDXs.size() ; i ++) pooe(i);
}

