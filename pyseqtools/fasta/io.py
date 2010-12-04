def seek_to_start(fh, max_iter=100):
    """Moves a file handle to the first ">" in a FASTA file.
    
    All lines that come before that are ignored.
    """
    if isinstance(fh, str):
        fh = open(fh, 'r')
    fh.seek(0)
    found = False
    for idx,line in enumerate(fh):
        if idx > max_iter:
            raise Exception("Seeking passed max_iter")
        line = line.strip()
        if line.startswith(">"):
            found = True
            break
    if not found:
        raise Exception("Seeking past end of file")
    fh.seek(0)
    while idx > 0:
        fh.next()
        idx = idx - 1
    return fh

