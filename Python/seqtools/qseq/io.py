from seqtools.qseq import QseqRead
from seqtools.io import xopen

"""
qseq format is tab delimted with the following columns:

  1. machine name (unique)
  2. run number (unique)
  3. lane number [1,8]
  4. tile number (positive integer)
  5. x coordinate of spot (can be negative)
  6. y coordinate of the spot (can be negative)
  7. index (positive number, should greater than zero, but most files are 0(?))
  8. read number (1 for single read, 2 for paired end)
  9. sequence (`.` instead of N)
  10. quality (phred) scores
  11. QC filter (1: passed QC, 0: no)
  
HISEQ2  195     8       1101    1402    1998    0       1       .AGAGCACACGTCTGA.C..A......CGATGTC...............CC     BOW]______[[[___BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB     0

"""
def parse(qseq, *args, **kwargs):
    fh = xopen(qseq, 'r')
    record = None
    for idx,line in enumerate(fh):
        qseq = QseqRead.from_record(line)
        
        