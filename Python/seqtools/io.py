import gzip

def xopen(x, mode='r', is_gzipped=None):
    """
    Replacement for the "open" function that can also open
    files that have been compressed with gzip. If the filename ends with .gz,
    the file is opened with gzip.open(). If it doesn't, the regular open()
    is used.
    closing -- whether to wrap the returned GzipFile with contextlib.closing
        to make it usable in 'with' statements. (TODO look for a nicer solution)
    
    Taken from cutadapt
    """
    if isinstance(x, file):
        if is_gzipped:
            return gzip.GzipFile(fileobj=x, mode='rb')
        else:
            return x
    
    if not isinstance(x, str):
        raise ValueError("Don't know how to handle objects of type: " + str(type(x)))
    
    if is_gzipped is None:
        is_gzipped = x.endswith('.gz')
    
    if is_gzipped:
        # if is_closing:
        #     return closing(gzip.open(filename, mode))
        # else:
        #     return gzip.open(filename, mode)
        if mode.find('b') == -1:
            mode = mode + 'b'
        return gzip.open(x, mode)
    else:
        return open(x, mode)
