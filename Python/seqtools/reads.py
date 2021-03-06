import re

from seqtools.utilities.enum import Enum
from seqtools.qualities import NGSQuality

NGSSequenceSpace = Enum('BASE', 'COLOR')

class NGSRead(object):
    """Base class for an object that represents a read from a sequencer"""
    
    def __init__(self, id, sequence=None, quality=None,
                 sequence_space=NGSSequenceSpace.BASE,
                 quality_type=NGSQuality.UNKNOWN):
        self.id = id
        self.sequence = sequence
        self.quality = quality
        self.sequence_space = sequence_space
        self.quality_type = quality_type
    
    def __len__(self):
        if self.sequence is None:
            l = 0
        else:
            l = len(self.sequence)
            if self.sequence_space == NGSSequenceSpace.COLOR:
                l = l - 1
        return l
    
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise TypeError("NGSRead objects only subscriptable by" \
                            "character %s" % str(key))
        return self.__dict__[key]
    
    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError("NGSRead objects only subscriptable by character")
        self.__dict__[key] = value
    
    def __repr__(self):
        """Basic representation of the read"""
        repr = "\n".join(["{'id' : '>%s'" % self.id,
                          " 'sequence' : '%s'" % self.sequence,
                          " 'quality' : '%s'}" % str(self.quality)])
        return repr
    
    def trim(self, n, side="left", minlength=4):
        """Number of basepairs to trim off of the given side."""
        if n <= 0:
            return
        self.sequence = trim_sequence(self.sequence, n, side, minlength,
                                      self.sequence_space)
        if self.quality is not None:
            self.quality = trim_quality(self.quality, n, side, minlength,
                                        self.quality_type)
        return self
    
###############################################################################
## Utilities
colorspace_homopolymer = re.compile(r"""^0+$""")
def trim_sequence(sequence, n, side='right', minlength=4,
                  sequence_space=NGSSequenceSpace.BASE):
    """Trim all types of sequences.
    
    TODO: Implement trimming non-homopolymers from Left in colorspace
    
    """
    if n <= 0:
        return sequence
    
    seqlen = len(sequence)
    if sequence_space == NGSSequenceSpace.COLOR:
        seqlen = seqlen - 1
    if seqlen - n < minlength:
        raise ValueError("Length would be less than %d" % minlength)
    
    ## Arugments check out: do the trimming
    if side == 'right':
        return sequence[:-n]
    
    ## Trimming from left is dangerous in colorspace
    if sequence_space == NGSSequenceSpace.BASE:
        return sequence[n:]

    if colorspace_homopolymer.match(sequence[1:1+n]):
        return sequence[0] + sequence[n+1:]
    else:
        raise NotImplementedError("Need to recode anchor base")

def trim_quality(quality, n, side='right', minlength=4,
                 quality_type=NGSQuality.UNKNOWN):
    ## Default SOLID Qualities are typically represented as a character string
    ## of numbers. Each quality is a whole number (ie. 21)
    if quality_type == NGSQuality.UNKNOWN:
        easy = isinstance(quality, list) or quality.find(' ') == -1
    else:
        easy = isinstance(quality, list) or quality_type != NGSQuality.SOLID
    
    if easy:
        qlen = len(quality)
        if qlen - n < minlength:
            raise ValueError("Length would be less than %d" % minlength)
        if side == 'right':
            return quality[:-n]
        else:
            return quality[n:]
    
    if side == 'right':
        while n > 0:
            quality = quality[:quality.rfind(' ')]
            n = n - 1
    else:
        while n > 0:
            quality = quality[quality.find(' ') + 1:]
            n = n - 1
    return quality
