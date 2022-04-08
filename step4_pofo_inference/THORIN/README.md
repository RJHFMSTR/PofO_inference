# Target Haplotype ORigin INference

## Overview




If you use the THORIN in your research work, please cite the following paper:

[Hofmeister et al., bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.11.03.467079v1)

## Documentaion

To run an example, use the following command line:

```
./bin/thorin_static -I test/related.chr20.vcf.gz -H test/unrelated.chr20.vcf.gz -G test/Trios.txt -M maps/chr20.b37.gmap.gz -R 20 -O output_test.txt
```

Options:

	**--help** outputs the list of possible options.
	```-I``` specifies the input .vcf file with all samples listed in the group file (-G).
	```-H``` specifies the input .vcf file with a set of fully unrelated samples.
	```-G``` specifies the group file.
	```-M``` specifies the genetic maps.
	```-O``` specifies the output file


The group file (```-G```) is formatted as follows:
* the first column specifies the target individual id
* the second column specifies the first group of relatives. It starts with ```G1=```, referring to ```Group 1```, followed by a ';'-separated list of relatives.
* the third column, if any, specifies the second group of relatives with ```G2=```.


The output file stores the copying probability of each group of relatives (```G1``` and ```G2``` if any) as well as the copying probability of unrelated individuals (```H```). It is formated as follows:
* columns 1-4 specifies the chromosome, position, index of the variant and centimorgan.
* each additional column specifies the probability of copying a target haplotype from a group of individuals.

* Ex: a column named ```NA20900_G1_0``` stores the probability of copying the haplotype 0 of the individuals NA20900 from the group G1.
* Ex: a column named ```NA20900_HOLE_1``` stores the probability of copying the haplotype 1 of the individuals NA20900 from the group H (unrelated individuals).

### Installation

A static version of THORIN ready to run is provided in ```./bin/thorin_static```




If you want to compile the code, make sure that the two following libraries are installed on your system:

	- HTSlib: A great C library for reading/writing high-throughput sequencing data.
	- BOOST: A free peer-reviewed portable C++ source libraries. SHAPEIT4 uses two specific BOOST libraries: iostreams and program_options.



Make sure that the following standard library flags can be used by g++ on your system:

	-lz, -lbz2 and -llzma for reading/writing compressed files.
	-lm for basic math operations.
	-lpthread for multi-threading.



Edit the makefile at lines 5-6 and 9-11 so that the following variables are correctly set up (look at the paths already there for an example):

    HTSLIB_INC (line 5): path to the HTSlib header files,
    HTSLIB_LIB (line 6): path to the static HTSlib library (file libhts.a),
    BOOST_INC (line 9): path to the BOOST header files (often /usr/include),
    BOOST_LIB_IO (line 10): path to the static BOOST iostreams library (file libboost_iostreams.a),
    BOOST_LIB_PO (line 11): path to the static BOOST program_options library (file libboost_iostreams.a),



Once all paths correctly set up, proceed with the THORIN compilation using ```make```.




## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

