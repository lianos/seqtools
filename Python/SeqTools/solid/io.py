from __future__ import absolute_import

from itertools import izip
import os,re

from SeqTools.fasta.io import seek_to_start

from SeqTools.solid.reads import SolidRead
from SeqTools.solid import dibase, convert

def parse(csfasta, csqual=None, sequence_convert=None,
          quality_convert=None):
    """Returns sequence and quality records one at a time via a generator.
    
    Arguments:
    - `csfasta`: The path to the csfasta file
    - `csqual`: The path to the quality file
    - `sequence_cnvert`: A list of function names used to convert the sequence
       of the reads. These functions are applied in the same order they are
       given.
    - `quality_convert` : A list of function names used to convert the quality
       scores of reads. These functions are applied in the same order they
       are given.
    
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
    
    ## Load the conversion functions
    if sequence_convert is not None and len(sequence_convert) > 0:
        if isinstance(sequence_convert, str):
            sequence_convert = [sequence_convert]
        sequence_convert = [seq_convert[wut] for wut in sequence_convert]
    
    if quality_convert is not None and len(quality_convert) > 0:
        if isinstance(quality_convert, str):
            quality_convert = [quality_convert]
        quality_convert = [qual_convert[wut] for wut in qual_convert]
    
    ## Determine if we are parsing qualities along with csfasta files or not
    if csqual is None:
        _parse = _parse_csfasta_only(csfasta)
    else:
        _parse = _parse_csfasta_with_qual(csfasta, csqual)
    
    ## Get down to business
    for record in _parse:
        if sequence_convert is not None:
            for sc in sequence_convert:
                record.sequence = sc(record.sequence)
        if quality_convert is not None:
            for qc in quality_convert:
                record.quality = qc(record.quality)
        yield record
    raise StopIteration("Parsing finished")

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

