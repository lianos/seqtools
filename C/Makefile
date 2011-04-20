# This assumes that libgtextutils, fastx-toolkit and their header files
# are in an "accessible" place, see:
# https://github.com/lianos/seqtools/wiki/FASTX-toolkit-Extensions
CXXFLAGS := -Wall -I. -L./sequtils -lsequtils -lgtextutils -lfastx
CXX := g++

FASTXX := fastxx_homopolymer_trimmer
SEQUTILS := sequtils/strings.cpp

all: fastxx
	@echo Making all

fastxx: fastxx_homopolymer_trimmer/fastxx_homopolymer_trimmer.cpp utils
	$(CXX) $(CXXFLAGS) fastxx_homopolymer_trimmer/fastxx_homopolymer_trimmer.cpp -o ../bin/fastxx_homopolymer_trimmer

utils: $(SEQUTILS)
	$(CXX) -c -Wall -I. sequtils/strings.cpp -o sequtils/sequtils.o
	ar rcs sequtils/libsequtils.a sequtils/sequtils.o

clean:
	-rm sequtils/*.a
	-rm sequtils/*.o

