#COMPILER MODE C++11
CXX=g++ -std=c++11


### Speficy your own path for HTSLIB and BOOST


#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
#HTSLIB_INC=$(HOME)/Tools/htslib-1.9
#HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
#BOOST_INC=/usr/include
#BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
#BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a



#COMPILER & LINKER FLAGS

#Best performance is achieved with this. Use it if running on the same plateform you're compiling, it's definitely worth it!
#CXXFLAG=-O3 -march=native

#Good performance and portable on most intel CPUs
CXXFLAG=-O3

#Portable version without avx2 (much slower)
#CXXFLAG=-O3

LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma

#SHAPEIT SOURCES & BINARY
BFILE=bin/pooer
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#COMPILATION RULES
all: $(BFILE)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
