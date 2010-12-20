from seqtools.solid import dibase

###############################################################################
## Sequence Conversions
def colorspace_to_basespace(x):
    """Convert color space to basespace"""
    return dibase.decodeSequence(x)

def basespace_to_colorspace(x):
    return dibase.encodeSequence(x)

###############################################################################
## Quality Conversions
##
## http://en.wikipedia.org/wiki/FASTQ_format
## 
##  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
##  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
##  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
##  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
##  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
##  |                         |    |        |                              |                     |
## 33                        59   64       73                            104                   126
## 
## S - Sanger        Phred+33,  raw reads typically (0, 40)
## X - Solexa        Solexa+64, raw reads typically (-5, 40)
## I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
## J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
##    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
##    (Note: See discussion above).

def quality_to_phred(quality_line, base=33, as_ascii=True):
    """Convert SOLiD quality to Sanger (base=33).
    
    Inspired from the galaxy::solid2sanger and cutadapt::quality_to_ascii
    functions.
    
    >>> quality_to_ascii("17 4 29 18")
    '2%>3'
    """
    phred = list()
    for x in quality_line.split(" "):
        x = int(x)
        if x < 0:
            x = 0
        phred.append(x + base)
    if as_ascii:
        phred = [chr(x) for x in phred]
        phred = ''.join(phred)
    return phred

def quality_to_sanger(quality, as_ascii=True):
    return quality_to_ascii(quality, 33, as_ascii=as_ascii)

def quality_to_integer(x):
    return [int(q) for q in x.split()]

def quality_to_string(x):
    return ' '.join([str(q) for q in x])

## Solexa->Sanger quality conversion table
## my @conv_table;
## for (-64..64) {
##   $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
## }
## ...
## $qual .= $conv_table[$_+64]; #intfq: SOLiD -> Standard

def quality_to_sanger(x):
    if isinstance(x, str):
        x = quality_to_integer(x)
    raise NotImplementedError()

def quality_to_illumina(x):
    x = quality_to_sanger(x)
    return [val - 33 for val in x]
