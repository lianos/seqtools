def sliding_window(sequence, n=3):
    """Returns a sliding window (of width n) over a character sequence"""
    total = len(sequence) - n + 1
    for i in range(total):
        yield sequence[i:i+n]
