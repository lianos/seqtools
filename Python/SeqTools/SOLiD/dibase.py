__doc__="""
    module for converting between 2-bp encoding and NT rep.
    
    Taken from Corona_Lite v.4.2.2
"""
import string

# change this if the probe labeling changes!
ENCODE_DICT = { \
    'AA':0,\
    'CC':0,\
    'GG':0,\
    'TT':0,\
    'AC':1,\
    'CA':1,\
    'GT':1,\
    'TG':1,\
    'AG':2,\
    'CT':2,\
    'GA':2,\
    'TC':2,\
    'AT':3,\
    'CG':3,\
    'GC':3,\
    'TA':3,\
    'A.':4,\
    'C.':4,\
    'G.':4,\
    'T.':4,\
    '.A':4,\
    '.C':4,\
    '.G':4,\
    '.T':4,\
    '.N':4,\
    'AN':4,\
    'CN':4,\
    'GN':4,\
    'TN':4,\
    'NA':4,\
    'NC':4,\
    'NG':4,\
    'NT':4,\
    'NN':4,\
    'N.':4,\
    '..':4\
}

# change this if the probe labeling changes!
DECODE_DICT = {
    'A0':'A',\
    'A1':'C',\
    'A2':'G',\
    'A3':'T',\
    'A4':'N',\
    'A.':'N',\
    'C0':'C',\
    'C1':'A',\
    'C2':'T',\
    'C3':'G',\
    'C4':'N',\
    'C.':'N',\
    'G0':'G',\
    'G1':'T',\
    'G2':'A',\
    'G3':'C',\
    'G4':'N',\
    'G.':'N',\
    'T0':'T',\
    'T1':'G',\
    'T2':'C',\
    'T3':'A',\
    'T4':'N',\
    'T.':'N',\
    'N0':'N',\
    'N1':'N',\
    'N2':'N',\
    'N3':'N',\
    'N.':'N'\
}

def decodeSequence( encodedSequence ):
    " Converts from 2-base encoding to sequence space "
    currentBase = encodedSequence[0]
    sequence = []
    for i in range(1,len(encodedSequence)):
        dibase = currentBase + encodedSequence[i]
        currentBase = DECODE_DICT[ dibase ]
        sequence.append( currentBase )
    return "".join(sequence)

def encodeSequence( sequence ):
    " Returns 2-base encoding for a nucleotide sequence "
    currentBase = sequence[0]
    codes = [ currentBase ]
    for i in xrange( 1, len(sequence) ):
        dibase = currentBase + sequence[i]
        codes.append( str(ENCODE_DICT[dibase]) )
        currentBase = sequence[i]
        
    return "".join(codes)
    
def doubleDecodeSequence( deSequence ):
    " Converts a double-encoded sequence to base space "
    table = string.maketrans( 'NACGT', '.0123' )
    eSequence = deSequence[0] + deSequence[1:].translate(table)
    return decodeSequence( eSequence )
    
def doubleEncodeSequence( eSequence ):
    " Converts a 2-base encoded sequence to a double encoded nucleotide sequence "
    table = string.maketrans( '.01234', 'NACGTN' )
    return eSequence.translate( table )
