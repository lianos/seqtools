from gzip import GzipFile
from seqtools.fastq import FastqRead
from seqtools.io import xopen

def parse(infile, *args, **kwargs):
    if isinstance(infile, file) or isinstance(infile, GzipFile):
        fh = infile
        do_close = False
    else:
        fh = xopen(infile, 'r')
        do_close = True
    record = None
    for idx,line in enumerate(fh):
        line = line.strip()
        step = idx % 4
        if step == 0:
            if record is not None:
                yield record
            record = FastqRead(id=line[1:])
        if step == 1:
            record.sequence = line
        elif step == 2:
            record.optional_id = line[1:]
        else:
            record.quality = line
    if do_close:
        fh.close()
