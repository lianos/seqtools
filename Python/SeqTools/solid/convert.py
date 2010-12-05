from SeqTools.solid import dibase

###############################################################################
## Sequence Conversions
def colorspace_to_basespace(x):
    """Convert color space to basespace"""
    return dibase.decodeSequence(x)

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
