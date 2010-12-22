import gzip

def xopen(filename, mode='r'):
    """
    Replacement for the "open" function that can also open
    files that have been compressed with gzip. If the filename ends with .gz,
    the file is opened with gzip.open(). If it doesn't, the regular open()
    is used.
    closing -- whether to wrap the returned GzipFile with contextlib.closing
        to make it usable in 'with' statements. (TODO look for a nicer solution)
    
    Taken from cutadapt
    """
    if filename.endswith('.gz'):
        # if is_closing:
        #     return closing(gzip.open(filename, mode))
        # else:
        #     return gzip.open(filename, mode)
        if mode.find('b') == -1:
            mode = mode + 'b'
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
