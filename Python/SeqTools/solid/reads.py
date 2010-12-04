from SeqTools import NGSRead, NGSQuality

class SolidRead(NGSRead):
    """Documentation for SolidRecord"""
    
    def __init__(self, id, sequence=None, quality=None,
                 sequence_space=NGSSequenceSpace.COLOR,
                 quality_type=NGSQuality.SOLID):
        super(SolidRead, self).__init__(id, sequence, quality, sequence_space,
                                        quality_type)
    
    def __repr__(self):
        """docstring for __repr__"""
        repr = "\n".join(["{'id' : '>%s'" % self.id,
                          " 'sequence' : '%s'" % self.sequence,
                          " 'quality' : '%s'}" % str(self.quality)])
        return repr
    
    def __len__(self):
        if self.in_colorspace and self.sequence is not None:
            return len(self.sequence) - 1
        else:
            return super(SolidRead,self).__len__()

###############################################################################
## Utilities
