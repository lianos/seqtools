LIB_FASTX = /Users/stavros/cBio/projects/fastx_toolkit/fastx/src/libfastx
LIB_GTEXTUTILS = /Users/stavros/cBio/projects/fastx_toolkit/gtextutils

FASTXX_SRC = fastxx_homopolymer_trimmer/fastxx_homopolymer_trimmer.cpp
CFLAGS = -I$(LIB_FASTX) -I$(LIB_GTEXTUTILS)
CC = g++ -Wall

fastxx: $(FASTXX_SRC)
	$(CC) $(CFLAGS) $(FASTXX_SRC)

all: fastxx
	@echo Making all



    