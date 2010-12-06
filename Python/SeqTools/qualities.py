import numpy as np
from SeqTools.utilities.enum import Enum

## Classes and code to primarily deal with gathering stats over NGS primary
## in a reasonably efficient (in memory) manner

NGSQuality = Enum('UNKNOWN', 'SOLEXA', 'ILLUMINA', 'SANGER', 'SOLID')

class QualityMatrix(dict):
    """Store quality score observations in a space efficient fashion.
    
    This object does not assume a specific range/order of the quality scores.
    The ``QualityMatrix.rows`` is used to keep track of which row a given
    quality score (represented as a string) is stored in.
    """
    
    def __init__(self, read_length=32, quality_type=NGSQuality.UNKNOWN,
                 encoding='ASCII'):
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
        """Record an observation of a single quality score at given position"""
        if position >= self.read_length:
            raise ValueError("Position of quality observation is out of bounds")
        if self.quality_type == NGSQuality.SOLID:
            quality = int(quality)
        self[quality][position] += 1
        self.nobs_at_position[position] += 1
    
    def observe_all(self, quality_sequence):
        """Observe all quality scores for a given record.
        
        ``quality_sequence`` should be enumerable. If it is string of
        space-separated integers (like SOLID quality scores from their fasta
        qual files), make sure to split() the string before passing it in here.
        
        """
        for idx,value in enumerate(quality_sequence):
            self.observe_at(idx, value)
    
    def as_matrix(self):
        order = sorted(self.keys())
        m = np.zeros((len(self), self.read_length))
        for idx,value in enumerate(order):
            m[idx,:] = self[value]
        return m
    
    def __str__(self):
        order = sorted(self.keys())
        output = list()
        for key in order:
            value = [str(key)]
            value.extend([str(x) for x in self[key]])
            output.append(' '.join(value))
        return '\n'.join(output)
    
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

