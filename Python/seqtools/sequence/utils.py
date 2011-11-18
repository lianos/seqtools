import string

__complement_table = {'dna' : string.maketrans('ACGTNacgtn', 'TGCANtgcan'),
                      'rna' : string.maketrans('ACGUNacgun', 'UGCANugcan')}

def complement(x, type='dna'):
    type = type.lower()
    if type not in ('dna', 'rna'):
        raise ValueError("reverse_complement only works with DNA or RNA")
    return x.translate(__complement_table[type])
    
def reverse_complement(x, type='dna'):
    return complement(x, type)[::-1]
