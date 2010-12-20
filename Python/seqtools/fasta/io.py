from seqtools.fasta import FastaRecord

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

def parse(fasta):
    """Parse a fasta file into a FastaRecord.
    
    You might consider using Biopython.
    """
    fh = seek_to_start(fasta)
    record = None
    for line in fh:
        line = line.strip()
        is_new = line.startswith(">")
        if is_new:
            if record is not None:
                record.value = "".join(record.value)
                yield record
            record = FastaRecord(line, list())
        else:
            record.value.append(line)
    fh.close()
    raise StopIteration("Parsing finished")
