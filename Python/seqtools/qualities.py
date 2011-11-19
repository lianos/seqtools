"""
.. module:: qualities
   :platform: Unix, Windows
   :synopsis: Tools to deal with (phred) quality scores fore reads

.. moduleauthor:: Steve Lianoglou <slianoglou@gmail.com>

How to handle PHRED/Solexa/Illumina scores
http://lists.open-bio.org/pipermail/biopython-dev/2009-February/005386.html

Phred offsets for common scores:
    Sanger              33
    Illumina 1.9+       33
    Illumina pre 1.9    64

"""

#import numpy as np
from seqtools.utilities.enum import Enum

## Classes and code to primarily deal with gathering stats over NGS primary
## in a reasonably efficient (in memory) manner

NGSQuality = Enum('UNKNOWN', 'SOLEXA', 'ILLUMINA', 'SANGER', 'SOLID')
NGSQualityEncoding = Enum('ASCII', 'INTEGER')

NGSQualityBaseOffset = {
    'illumina' : 33,
    'sanger' : 33,
    'solexa' : 64,
    'illumina_old' : 64
}

def match_phred_offset(x):
    val = x
    if isinstance(x, str):
        try:
            val = int(x)
        except ValueError:
            try:
                val = NGSQualityBaseOffset[x.lower()]
            except KeyError:
                raise ValueError("Unknown phred base offset: " + str(x))
    if not isinstance(val, int):
        raise ValueError("Unknown input type for match " + str(type(val)))
    if val not in NGSQualityBaseOffset.values():
        raise ValueError("Unknown phred base offset: " + str(val))
    return val

## Reference for converting qualities:
## http://jumpgate.caltech.edu/wiki/QSeq
## Qualities are assumed to be ASCII-encoded as chr(qual + base)
try:
    import numpy
    def convert_quality(qual, inbase=NGSQualityBaseOffset['solexa'],
                        outbase=NGSQualityBaseOffset['sanger']):
        if inbase == outbase:
            return qual
        quality = numpy.asarray(qual, 'c')
        quality.dtype = numpy.uint8
        quality -= inbase - outbase
        quality.dtype = '|S1'
        return quality.tostring()

except ImportError:
    def convert_quality(qual, inbase=NGSQualityBaseOffset['solexa'],
                        outbase=NGSQualityBaseOffset['sanger']):
        if inbase == outbase:
            return qual
        quality = (chr(ord(x) - (inbase - outbase)) for x in qual)
        return ''.join(quality)

class QualityMatrix(dict):
    """Store quality score observations in a space efficient fashion.
    
    This object does not assume a specific range/order of the quality scores.
    """
    
    def __init__(self, read_length=32, quality_type=NGSQuality.UNKNOWN,
                 encoding='NGSQualityEncoding.ASCII'):
        """Build a QualityMatrix for the given ``read_length``.
        
        Args:
            read_length (int): The (maximum) length of the reads in the 
                experiment
            quality_type (int): The ``NGSQuality`` ``Enum`` specifying the
                "space" the quality scores are in (``UNKNOWN``, ``SOLEXA``,
                etc.).
            encoding (int): The ``NGSQualityEncoding`` Enum specifying whether
                the quality score is ascii or integer.
        
        """
        super(QualityMatrix, self).__init__()
        self.read_length = read_length
        self.quality_type = quality_type
        self.nobs_at_position = [0] * read_length
        self.encoding = encoding
    
    def __missing__(self, key):
        value = [0] * self.read_length
        self[key] = value
        return value
    
    def observe_at(self, position, quality):
        """Record an observation of a single quality score at given position.
        
        Args:
            position (int): The position (0-based index) of the read the
                quality is observed at.
            quality : The quality score observed at the position.
        
        Returns:
            None
        """
        if position >= self.read_length:
            raise ValueError("Position of quality observation is out of bounds")
        if self.quality_type == NGSQuality.SOLID:
            quality = int(quality)
        self[quality][position] += 1
        self.nobs_at_position[position] += 1
    
    def observe_all(self, quality_sequence):
        """Observe all quality scores for a given record.
        
        Args:
            quality_sequence (iterable): If it is string of space-separated
            integers (like SOLID quality scores from their fasta qual files),
            make sure to ``split()`` the string before passing it in here.
            Otherwise a "standard" FASTQ/Phred score can be passed, eg.
            ``aR`\HZFW[KT^GUVF[STGFQ\[Q___TY]_T_BB`` would be kosher
        
        Returns:
            None
        """
        for idx,value in enumerate(quality_sequence):
            self.observe_at(idx, value)
    
    def as_matrix(self):
        """Converts the ``QualityMatrix`` to a NumPy matrix.
        
        Returns:
            A ``NumPy`` matrix of the quality scores.
        """
        order = self.matrix_keys()
        m = np.zeros((len(self), self.read_length))
        for idx,value in enumerate(order):
            m[idx,:] = self[value]
        return m
    
    def matrix_keys(self):
        """Defines the order that the rows of matrix are printed"""
        return sorted(self.keys())
    
    def __str__(self):
        order = sorted(self.keys())
        output = list()
        for key in order:
            value = [str(key)]
            value.extend([str(x) for x in self[key]])
            output.append(' '.join(value))
        return '\n'.join(output)
    
    ###########################################################################
    ## Bit rot
    def __initialize_qualities(self):
        """Initializes the dict to the known quality scores.
        
        TODO: Implement __initialize_qualities
        """
        ranges = {
          NGSQuality.UNKNOWN : """""",
          ## 59 - 104
          NGSQuality.SOLEXA : """;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""",
          ## 64 - 104
          NGSQuality.ILLUMINA : """@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh""",
          NGSQuality.SANGER : """!"#$%&'()*+,-./0123456789:;""", ## 33 - 73
          NGSQuality.SOLID : """""",
        }
        ranges = ranges[self.quality_type]
        if self.encoding != 'ASCII':
            ## TODO: Deal with char <-> Int conversions
            pass
        [self[str(x)] for x in ranges]

