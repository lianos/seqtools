# This assumes that libgtextutils, fastx-toolkit and their header files
# are in an "accessible" place, see:
# https://github.com/lianos/seqtools/wiki/FASTX-toolkit-Extensions
CFLAGS = -Wall -I. -L./sequtils -lsequtils -lgtextutils -lfastx
CC = g++

all: fastxx
	@echo Making all

fastxx: fastxx_homopolymer_trimmer/fastxx_homopolymer_trimmer.cpp utils
	$(CC) $(CFLAGS) fastxx_homopolymer_trimmer/fastxx_homopolymer_trimmer.cpp -o ../bin/fastxx_homopolymer_trimmer

utils: sequtils/strings.cpp
	$(CC) -c -Wall -I. sequtils/strings.cpp -o sequtils/sequtils.o
	ar rcs sequtils/libsequtils.a sequtils/sequtils.o

clean:
	-rm sequtils/*.a
	-rm sequtils/*.o

