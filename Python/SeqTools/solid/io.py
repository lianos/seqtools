from __future__ import absolute_import

from itertools import izip
import os,re

from ..fasta.io import seek_to_start
from . import dibase
from . import convert

class SolidRead(object):
    """Documentation for SolidRecord"""
    
    def __init__(self, _id, sequence=None, quality=None, quality_type='original'):
        self.id = _id
        self.sequence = sequence
        self.quality = quality
        self.quality_type = quality_type
        self.in_colorspace = True
    
    def __trim_quality(self, n, side='left'):
        if side == 'left':
            if self.quality_type == 'original':
                while n > 0:
                    self.quality = self.quality[(self.quality.find(' ') + 1):]
                    n = n - 1
            elif self.quality_type == 'integer':
                self.quality = self.quality[n:]
        else:
            raise Exception("Trimming from right is easy")
    
    def trim(self, n, side="left"):
        if side == 'left':
            self.sequence = self.sequence[n:]
        else:
            self.sequence = self.sequence[:-n]
        self.__trim_quality(n, side)
    
    def __repr__(self):
        """docstring for __repr__"""
        repr = "\n".join(["{'id' : '>%s'" % self.id,
                          " 'sequence' : '%s'" % self.sequence,
                          " 'quality' : '%s'}" % str(self.quality)])
        return repr

# END : Class SolidRecord

def _parse_csfasta_only(csfasta):
    ffh = seek_to_start(csfasta)
    for fasta in ffh:
        fasta = fasta.strip()
        fnew = fasta.startswith(">")
        if fnew:
            fid = fasta[1:]
            # record = {'id' : fid, 'sequence' : None, 'quality' : None}
            record = SolidRead(fid)
        else:
            record.sequence = fasta
            yield record
    ffh.close()
    raise StopIteration("Parsing finished")

def _parse_csfasta_with_qual(csfasta, qual):
    ffh = seek_to_start(csfasta)
    qfh = seek_to_start(qual)
    
    for fasta,qual in izip(ffh, qfh):
        fasta = fasta.strip()
        fnew = fasta.startswith(">")
        qual = qual.strip()
        qnew = qual.startswith(">")
        
        if fnew ^ qnew:  # xor
            raise Exception("Mismatched files: id and info not aligned")
        elif fnew and qnew:
            fid = fasta[1:]
            qid = qual[1:]
            if fid != qid:
                raise Exception("Mismatched files: id's do not match")
            # record = {'id' : fid, 'sequence' : None, 'quality' : None}
            record = SolidRead(fid)
        else:
            record.sequence = fasta
            record.quality = qual
            yield record
    ffh.close()
    qfh.close()
    raise StopIteration("Parsing finished")

def parse(csfasta, csqual=None, sequence_convert=None,
          quality_convert=None):
    """Returns sequence and quality records one at a time via a generator.
    
    Arguments:
    - `csfasta`: The path to the csfasta file
    - `csqual`: The path to the quality file
    - `comment_char`: Skip lines with these characters. Assumes they only appear
    at the start of the file
    
    """
    seq_convert = {
        'basespace' : convert.colorspace_to_basespace,
    }
    qual_convert = {
        'integer'   : convert.quality_to_integer,
        'sanger'    : convert.quality_to_sanger,
        'solexa'    : convert.quality_to_sanger,
        'illumina'  : convert.quality_to_illumina,
    }
    
    if sequence_convert is not None:
        if isinstance(sequence_convert, str):
            sequence_convert = [sequence_convert]
        sequence_convert = [seq_convert[wut] for wut in sequence_convert]
    if quality_convert is not None:
        if isinstance(quality_convert, str):
            quality_convert = [quality_convert]
        quality_convert = [qual_convert[wut] for wut in qual_convert]
    
    if csqual is None:
        _parse = _parse_csfasta_only(csfasta)
    else:
        _parse = _parse_csfasta_with_qual(csfasta, csqual)
    
    for record in _parse:
        if sequence_convert is not None:
            for sc in sequence_convert:
                record.sequence = sc(record.sequence)
        if quality_convert is not None:
            for qc in quality_convert:
                record.quality = qc(record.quality)
        yield record
    raise StopIteration("Parsing finished")



## strip N quality scores
def strip_qual(qual, n, side='left'):
    if side == 'left':
        if isinstance(qual, str):
            while n > 0:
                qual = qual[(qual.find(' ') + 1):]
                n = n -1
        else:
            qual = qual[n:]
    else:
        raise Exception("Trimming from right is easy")
    return qual


_all_bases_colorspace = re.compile("^[ACGT]\d+$")
_all_bases_basespace = re.compile("^[ACGT]+$")
def filter_all_bases(record, in_basespace=False):
    if in_basespace:
        regex = _all_bases_basespace
    else:
        regex = _all_bases_colorspace
    return regex.match(record['sequence']) is not None


def filter_cleanseq(csfasta, csqual, keep=re.compile("^T\d+$"),
                    in_basespace=False):
    fname = '.'.join(['filter', os.path.basename(csfasta)])
    qname = '.'.join(['filter', os.path.basename(csqual)])
    fout = open(os.path.join(os.path.dirname(csfasta), fname), 'w')
    qout = open(os.path.join(os.path.dirname(csqual), qname), 'w')
    for idx,record in enumerate(parse(csfasta, csqual, convert_quality=False)):
        if keep.match(record['sequence']) is not None:
            fout.write(">" + record['id'] + "\n")
            fout.write(record['sequence'] + "\n")
            qout.write(">" + record['id'] + "\n")
            qout.write(record['quality'] + "\n")
    fout.close()
    qout.close()


def stripTrun(csfasta, csqual):
    fname = '.'.join(['stripT', os.path.basename(csfasta)])
    qname = '.'.join(['stripT', os.path.basename(csqual)])
    fout = open(os.path.join(os.path.dirname(csfasta), fname), 'w')
    qout = open(os.path.join(os.path.dirname(csqual), qname), 'w')
    regex = re.compile("T.*?[^0]")
    
    for idx,record in enumerate(parse(csfasta, csqual, convert_quality=False)):
        m = regex.match(record['sequence'])
        if m is not None:
            start = m.end() - 1
            fout.write(">" + record['id'] + "\n")
            fout.write('T' + record['sequence'][start:] + "\n")
            
            qout.write(">" + record['id'] + "\n")
            # qual = record['quality'][start:]
            # qout.write(' '.join([str(x) for x in qual]) + "\n")
            qout.write(strip_qual(record['quality'], start) + "\n")
    fout.close()
    qout.close()
