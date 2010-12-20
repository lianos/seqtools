import copy

from seqtools import NGSRead, NGSSequenceSpace, NGSQuality
from seqtools.solid import convert

class SolidRead(NGSRead):
    """Documentation for SolidRecord"""
    
    def __init__(self, id, sequence=None, quality=None,
                 sequence_space=NGSSequenceSpace.COLOR,
                 quality_type=NGSQuality.SOLID):
        super(SolidRead, self).__init__(id, sequence, quality, sequence_space,
                                        quality_type)
    
    def to_basespace(self):
        obj = copy.deepcopy(self)
        if obj.sequence_space == NGSSequenceSpace.COLOR:
            obj.sequence = convert.colorspace_to_basespace(obj.sequence)
        return obj

###############################################################################
## Utilities
