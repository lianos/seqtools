import copy

from seqtools import NGSRead, NGSSequenceSpace, NGSQuality

"""
qseq format is tab delimted with the following columns:

  1. machine name (unique)
  2. run number (unique)
  3. lane number [1,8]
  4. tile number (positive integer)
  5. x coordinate of spot (can be negative)
  6. y coordinate of the spot (can be negative)
  7. index (positive number, should greater than zero, but most files are 0(?))
  8. read number (1 for single read, 2 for paired end)
  9. sequence (`.` instead of N)
  10. quality (phred) scores
  11. QC filter (1: passed QC, 0: no)
"""
class QseqRead(NGSRead):
    
    def __init__(self, machine="", run=-1, lane=-1, tile=-1, x=-1,
                 y=-1, index=0, read_no=1, sequence=None, quality=None,
                 qc=0, quality_type=NGSQuality.SANGER):
        self.machine = machine
        self.run = int(run)
        self.lane = int(lane)
        self.tile = int(tile)
        self.x = int(x)
        self.y = int(y)
        self.index = int(index)
        self.read_no = int(read_no)
        self.sequence = sequence
        self.sequence_space=NGSSequenceSpace.BASE
        self.quality = quality
        self.quality_type=quality_type
        self.qc = bool(qc)
        self.id = self.fastq_id()
    
    @staticmethod
    def from_record(line, sep="\t", convert_ambiguous=True):
        """Generate a QseqRead from a line in a (normal) qseq file
        
        qseq format is tab delimted with the following columns:

          1. machine name (unique)
          2. run number (unique)
          3. lane number [1,8]
          4. tile number (positive integer)
          5. x coordinate of spot (can be negative)
          6. y coordinate of the spot (can be negative)
          7. index (positive number, should greater than zero, but most files are 0(?))
          8. read number (1 for single read, 2 for paired end)
          9. sequence (`.` instead of N)
          10. quality (phred) scores
          11. QC filter (1: passed QC, 0: no)
        
        """
        info = line.strip().split(sep)
        read = QseqRead(machine=info[0], run=info[1], lane=info[2], tile=info[3],
                        x=int(info[4]), y=info[5], index=info[6], read_no=info[7],
                        sequence=info[8], quality=info[9], qc=info[10])
        return read

    def fastq_id(self):
        """Generate a fastq id for this record
        
        A fastq id looks like:
          HISEQ2:171:D0DDRABXX:8:1107:10531:157270
          machine:run:???:lane:tile:x:y
        """
        return "%s:%d:%s:%d:%d:%d:%d" % (self.machine, self.run, '???',
                                         self.lane, self.tile, self.x, self.y)

    def to_fastq(self):
        """Convert the qseq read to fastq
        """
        fqid = self.fastq_id()
        return FastqRead(fqid, self.sequence, quality=self.quality,
                         sequence_space=self.sequence_space,
                         quality_type=self.quality_type)
