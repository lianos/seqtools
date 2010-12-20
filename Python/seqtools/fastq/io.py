from seqtools.fastq import FastqRead

def parse(fastq, *args, **kwargs):
    fh = open(fastq, 'r')
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
    fh.close()
