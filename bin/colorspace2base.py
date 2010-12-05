#!/usr/bin/env python

from optparse import OptionParser
import os,sys

from SeqTools.solid import SolidRead
from SeqTools.solid.convert import colorspace_to_basespace
from SeqTools.solid.io import parse

if __name__ == '__main__':
    usage = "usage: %prog [options] CSFASTA\n\n" \
            "Converts a fasta file from colorspace to base space.\n" \
            "The result is sent to STDOUT"
    parser = OptionParser(usage=usage)
    (options,args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Path to csfasta file required")
    csfasta = args[0]
    
    if not os.path.isfile(csfasta):
        parser.error("FASTA file can not be read")
    
    for idx,read in enumerate(parse(csfasta)):
        # read = read.to_basespace()
        sys.stdout.write(">%s\n" % read.id)
        sys.stdout.write("%s\n" % colorspace_to_basespace(read.sequence))
    
    