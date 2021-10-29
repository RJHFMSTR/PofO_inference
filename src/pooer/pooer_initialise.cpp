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

#include <io/genotype_reader.h>
#include <io/gmap_reader.h>

void pooer::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step1: Read input files
	genotype_reader readerG(H, V, options["region"].as < string > ());
	if (options.count("hole")) {
		readerG.scanGenotypes(options["input"].as < string > (), options["hole"].as < string > ());
		readerG.allocateGenotypes();
		readerG.readGenotypes(options["input"].as < string > (), options["hole"].as < string > ());
	} else {
		readerG.scanGenotypes(options["input"].as < string > ());
		readerG.allocateGenotypes();
		readerG.readGenotypes(options["input"].as < string > ());
	}

	//step2: read grouping
	H.readGrouping(options["group"].as < string > ());

	//step3: Read and initialise genetic map
	gmap_reader readerGM;
	readerGM.readGeneticMapFile(options["map"].as < string > ());
	V.setGeneticMap(readerGM);
	M.initialise(V, options["effective-size"].as < int > (), readerG.n_samples*2);
}
