from collections import Counter
from seqtools.utilities import sliding_window

def frequency(sequence, k=3):
    freq = Counter()
    for kmer in sliding_window(sequence, k):
        freq[kmer] += 1
    return freq
