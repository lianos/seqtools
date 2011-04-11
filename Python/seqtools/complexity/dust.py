"""Dust scores are a semi-standard way to score complexity of sequence

See breakdown of dust sore here:
https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2009-February/000170.html
"""

import seqtools.sequence.kmer as kmer

def stats(sequence, k=3, square=True):
    """Returns the complexity score along with the number of unique kmers"""
    kfreq = kmer.frequency(sequence, k)
    if square:
        score = sum([(x - 1)**2 for x in kfreq.values()])
    else:
        score = sum([(x - 1) for x in kfreq.values()])
    return (score, len(kfreq))

def score(sequence, k=3, square=True):
    stats(sequence, k, squre)[0]

def max_score(sequence, k=3, square=True):
    ans = len(sequence) - k
    if square:
        ans = ans ** 2
    return ans

