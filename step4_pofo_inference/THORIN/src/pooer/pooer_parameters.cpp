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

void pooer::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread,T", bpo::value<int>()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value< string >(), "Genotypes to be phased in VCF/BCF format")
			("region,R", bpo::value< string >(), "Target region");
	
	bpo::options_description opt_group ("Group related files");
	opt_input.add_options()
			("group,G", bpo::value< string >(), "Group assignement for target samples")
			("hole,H", bpo::value< string >(), "Common group for all target samples");

	bpo::options_description opt_hmm ("HMM parameters");
	opt_hmm.add_options()
			("map,M", bpo::value< string >(), "Genetic map")
			("effective-size", bpo::value<int>()->default_value(15000), "Effective size of the population");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "HMM sumation per group")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_group).add(opt_hmm).add(opt_output);
}

void pooer::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("POOER");
	vrb.bullet("Author        : Olivier DELANEAU & Robin Hofmeister, University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com, robin@unil.ch");
	vrb.bullet("Version       : 1.1");
	vrb.bullet("Run date      : " + tac.date());
}

void pooer::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region or chromosome to phase using --region");

	if (!options.count("map"))
		vrb.error("You must specify a genetic map file using --map");

	if (!options.count("output"))
		vrb.error("You must specify an output file with --output");

	if (!options.count("group"))
		vrb.error("You must specify one input group file using --group");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");

	if (!options["thread"].defaulted() && !options["seed"].defaulted())
		vrb.warning("Using multi-threading prevents reproducing a run by specifying --seed");

	if (!options["effective-size"].defaulted() && options["effective-size"].as < int > () < 1)
		vrb.error("You must specify a positive effective size");
}

void pooer::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	vrb.bullet("Hole VCF      : [" + options["hole"].as < string > () + "]");
	vrb.bullet("Grouping      : [" + options["group"].as < string > () + "]");
	vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	vrb.bullet("Output        : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void pooer::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	vrb.bullet("HMM     : Recombination rates given by genetic map");
}
