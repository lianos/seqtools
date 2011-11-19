from seqtools.io import xopen

"""
Bowtie default output is almost sam like:
http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output

bowtie outputs one alignment per line. Each line is a collection of 8 fields separated by tabs; from left to right, the fields are:

1. Name of read that aligned

2. Reference strand aligned to, + for forward strand, - for reverse

3. Name of reference sequence where alignment occurs, or numeric ID if no name
   was provided

4. 0-based offset into the forward reference strand where leftmost character of
   the alignment occurs

5. Read sequence (reverse-complemented if orientation is -).

   If the read was in colorspace, then the sequence shown in this column is the
   sequence of decoded nucleotides, not the original colors. See the Colorspace
   alignment section for details about decoding. To display colors instead, use
   the --col-cseq option.

6. ASCII-encoded read qualities (reversed if orientation is -). The encoded
   quality values are on the Phred scale and the encoding is ASCII-offset by 33
   (ASCII char !).

   If the read was in colorspace, then the qualities shown in this column are the
   decoded qualities, not the original qualities. See the Colorspace alignment
   section for details about decoding. To display colors instead, use the
   --col-cqual option.

7. If -M was specified and the prescribed ceiling was exceeded for this read,
   this column contains the value of the ceiling, indicating that at least that
   many valid alignments were found in addition to the one reported.

   Otherwise, this column contains the number of other instances where the same
   sequence aligned against the same reference characters as were aligned
   against in the reported alignment. This is not the number of other places
   the read aligns with the same number of mismatches. The number in this
   column is generally not a good proxy for that number (e.g., the number in
   this column may be '0' while the number of other alignments with the same
   number of mismatches might be large).

8. Comma-separated list of mismatch descriptors. If there are no mismatches in
   the alignment, this field is empty. A single descriptor has the format
   offset:reference-base>read-base. The offset is expressed as a 0-based offset

"""
from seqtools import NGSRead, NGSSequenceSpace, NGSQuality

def parse(infile, *args, **kwargs):
    fh = xopen(infile, 'r')
    record = None
    for idx,line in enumerate(fh):
        yield BowtieRead.from_record(line)
    fh.close()

class BowtieRead(NGSRead):
    
    def __init__(self, id, sequence, quality=None, strand='*', seqname="*", start=0,
                 ceiling=0, mismatch=""):
        self.id = id
        self.sequence = sequence
        self.quality = quality
        self.strand = strand
        self.seqname = seqname
        self.start = int(start)
        self.ceiling = int(ceiling)
        self.mismatch = mismatch
    
    @staticmethod
    def from_record(line, sep="\t", convert_ambiguous=True):
        line = line.strip().split(sep)
        read_name = line[0]
        strand = line[1]
        seqname = line[2]
        start = line[3]
        seq = line[4]
        qual = line[5]
        ceiling = line[6]
        mismatch = "" if len(line) < 8 else line[7]
        return BowtieRead(read_name, seq, qual, strand, seqname,
                          start, ceiling, mismatch)

    def __str__(self):
        return "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\n" % \
          (self.id, self.strand, self.seqname, self.start, self.sequence,
           self.quality, self.ceiling, self.mismatch)
