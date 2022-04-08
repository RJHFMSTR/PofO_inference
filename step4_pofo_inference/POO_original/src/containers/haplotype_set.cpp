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
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_site = 0;
	n_hap = 0;
}

void haplotype_set::readGrouping(string ifile) {
	vrb.title("Reading target groups");
	string buffer;
	vector < string > tokens, tokensG, tokensS;
	map < string, int > :: iterator it;
	input_file fd(ifile);
	if (fd.fail()) vrb.error("Cannot open [" + ifile + "]");
	//getline(fd, buffer, '\n');
	while (getline(fd, buffer, '\n')) {
		if (stb.split(buffer, tokens) >= 2) {
			it = IDs_map.find(tokens[0]);
			if (it == IDs_map.end()) vrb.bullet("Sample [" + tokens[0] + "] unfound !");
			else {
				targetIDXs.push_back(it->second);
				sourceIDXs.push_back(vector < vector < int > > ());
				groupIDs.push_back(vector < string > ());

				for (int g = 1 ; g < tokens.size() ; g++) {
					stb.split(tokens[g], tokensG, "=");
					if (tokensG.size() != 2) vrb.error("Could not extract a group ID for [" + tokens[0] + "]");
					groupIDs.back().push_back(tokensG[0]);
					stb.split(tokensG[1], tokensS, ";");
					sourceIDXs.back().push_back(vector < int > ());
					for (int t = 0 ; t < tokensS.size() ; t ++) {
						it = IDs_map.find(tokensS[t]);
						if (it == IDs_map.end()) vrb.bullet("Sample [" + tokens[0] + "] has an undefined source sample [" + tokensS[t] + "] being unfound !");
						else sourceIDXs.back().back().push_back(it->second);
					}
				}
			}
		}
	}
	fd.close();

	hasHole = (IDs_map.size()*2 != n_hap);
	copyingProbabilities = vector < vector < float > > (2*targetIDXs.size());
	for (int t =0 ; t < targetIDXs.size() ; t ++) {
		copyingProbabilities[2*t+0] = vector < float > ((hasHole + sourceIDXs[t].size())*n_site*2, 0.0);
		copyingProbabilities[2*t+1] = vector < float > ((hasHole + sourceIDXs[t].size())*n_site*2, 0.0);
	}
}

void haplotype_set::writeCopyingProbabilities(string ofile, variant_map & V) {
	output_file fd (ofile);
	vrb.title("Writing Copying Probabilities in[" + ofile + "]");

	//Write file header
	fd << "#CHROM\tPOS\tIDX\tCM";
	for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
		//Hap0
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) fd << "\t" << IDs[targetIDXs[t]] << "_" << groupIDs[t][g] << "_0";
		if (hasHole) fd << "\t" << IDs[targetIDXs[t]] << "_HOLE_0";
		//Hap1
		for (int g = 0 ; g < sourceIDXs[t].size() ; g ++) fd << "\t" << IDs[targetIDXs[t]] << "_" << groupIDs[t][g] << "_1";
		if (hasHole) fd << "\t" << IDs[targetIDXs[t]] << "_HOLE_1";
		
	}
	fd << endl;

	//Write file body
	for (int l = 0 ; l  < n_site ; l ++) {
		fd << V.vec_pos[l]->chr << "\t" << V.vec_pos[l]->bp << "\t" << l << "\t" << stb.str(V.vec_pos[l]->cm, 3);
		for (int t = 0 ; t < sourceIDXs.size() ; t ++) {
			int ngroups = sourceIDXs[t].size() + hasHole;
			for (int g = 0 ; g < ngroups ; g ++) fd << "\t" << stb.str(copyingProbabilities[2*t+0][ngroups*l+g], 2);
			for (int g = 0 ; g < ngroups ; g ++) fd << "\t" << stb.str(copyingProbabilities[2*t+1][ngroups*l+g], 2);
		}
		fd << endl;
	}
	fd.close();
}

