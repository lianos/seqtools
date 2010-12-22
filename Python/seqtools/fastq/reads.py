import copy

from seqtools import NGSRead, NGSSequenceSpace, NGSQuality

class FastqRead(NGSRead):
    """Documentation for SolidRecord"""
    
    def __init__(self, id, sequence=None, optional_id=None, quality=None,
                 sequence_space=NGSSequenceSpace.BASE,
                 quality_type=NGSQuality.SANGER):
        super(FastqRead, self).__init__(id, sequence, quality, sequence_space,
                                        quality_type)
        self.optional_id=optional_id
    
    def to_basespace(self):
        obj = copy.deepcopy(self)
        if obj.sequence_space == NGSSequenceSpace.COLOR:
            obj.sequence = convert.colorspace_to_basespace(obj.sequence)
        return obj
    
    def __str__(self):
        return "@%s\n%s\n+\n%s\n" % (self.id, self.sequence, self.quality)

###############################################################################
## Utilities
