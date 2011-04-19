# INC_FASTX = /Users/stavros/cBio/projects/fastx_toolkit/fastx/src/libfastx
# INC_GTEXTUTILS = /Users/stavros/cBio/projects/fastx_toolkit/gtextutils
# LIB_GTEXTUTILS = /usr/local

INC_FASTX = /Volumes/Sequel/tools/nextgen_sequencing/FASTX-toolkit/fastx_toolkit-0.0.13
INC_LIBGTEXTUTILS = /Volumes/Sequel/tools/nextgen_sequencing/FASTX-toolkit/fastx_toolkit-0.0.13
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

    