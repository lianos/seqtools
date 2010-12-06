import numpy as np
from SeqTools.utilities.enum import Enum

NGSQuality = Enum('UNKNOWN', 'SOLEXA', 'ILLUMINA', 'SANGER', 'SOLID')

class QualityScoreIndex(dict):
    """Helper class used to return the row a given quality score appears"""
    
    def __init__(self, quality_type=NGSQuality.UNKNOWN, as_int=None):
        super(QualityScoreIndex, self).__init__()
        self.quality_type=quality_type
        if as_int is None:
            if self.quality_type in [NGSQuality.SOLID]:
                as_int = True
            else:
                as_int = False
        self.as_int = as_int
        self.__initialize_qualities()
        
    def __missing__(self, key):
        if self.as_int:
            key = int(key)
        value = len(self)
        self[key] = value
        return value
    
    def __initialize_qualities(self):
        """Initializes the dict to the known quality scores"""
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
        if self.as_int:
            ## TODO: Deal with char <-> Int conversions
            pass
        [self[str(x)] for x in ranges]
    

class QualityMatrix(object):
    """Store quality score observations in a space efficient fashion.
    
    This object does not assume a specific range/order of the quality scores.
    The ``QualityMatrix.rows`` is used to keep track of which row a given
    quality score (represented as a string) is stored in.
    """
    
    def __init__(self, read_length=32, quality_type=NGSQuality.UNKNOWN):
        self.read_length = read_length
        self.rows = QualityScoreIndex(quality_type)
        self.qmatrix = list()
        self.nobs_at_position = [0] * read_length
    
    def quality_distribution(self, quality):
        """Returns row representing frequency of quality score at each BP"""
        qidx = self.rows[quality]
        if qidx >= len(self.qmatrix):
            self.qmatrix.append([0] * self.read_length)
        return self.qmatrix[qidx]
    
    def observe_at(self, position, quality):
        """Record an observation of a quality score at given position"""
        if position >= self.read_length:
            raise ValueError("Position of quality observation is out of bounds")
        d = self.quality_distribution(quality)
        d[position] += 1
        self.nobs_at_position[position] += 1
    
    def observe_all(self, quality_sequence):
        for idx,value in enumerate(quality_sequence):
            self.observe_at(idx, value)
        
    def __str__(self):
        order = self.rows.keys()
        if self.rows.as_int:
            order = sorted(order)
        output = list()
        for value in order:
            row = self.quality_distribution(value)
            value = [str(value)]
            value.extend([str(x) for x in row])
            output.append(' '.join(value))
        return '\n'.join(output)
