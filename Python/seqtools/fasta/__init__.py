import re

_id_regex = re.compile(r"""\S+\s+""")
class FastaRecord(object):
    """Represent a FASTA object.
    
    A record like so ... ::
    
        >something goes here
        BLAHBLAHBLABH
        
    ... is represented by a FastaRecord ``fr`` like so::
    
        fr.id = "something"
        fr.description = "goes here"
        fr.value = "BLAHBLAHBLABH"
    
    """
    
    @staticmethod
    def parse_id(id):
        """Parse id string into its id and description"""
        id = id.strip()
        if id.startswith(">"):
            id = id[1:]
        if len(id) == 0:
            raise ValueError("Not a valid id string")
        m = _id_regex.match(id)
        if m is None:
            d = {"id" : id, "description" : ""}
        else:
            end = m.end()
            d = {"id" : id[:end-1].strip(), "description" : id[end:]}
        return d
        
    def __init__(self, id, sequence="", description="", *args, **kwargs):
        _id = FastaRecord.parse_id(id)
        if len(description) == 0:
            description = _id['description']
        self.id = _id['id']
        self.sequence = sequence
        self.description = description
    
    def __repr__(self):
        """Basic representation of the FASTSA Record"""
        repr = "\n".join(["{'id' : '>%s %s'" % (self.id, self.description),
                          self.sequence])
        return repr
    
    def __str__(self):
        if len(self.description):
            return ">%s %s\n%s\n" % (self.id, self.description, self.sequence)
        else:
            return ">%s\n%s\n" % (self.id, self.sequence)

#
from seqtools.fasta.io import parse

    
# END : Class FastaRecord
