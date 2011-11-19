"""Tools used to deal with raw NGS sequence files.

.. moduleauthor:: Steve Lianoglou <slianoglou@gmail.com>

"""
from seqtools.reads import NGSRead, NGSSequenceSpace
from seqtools.qualities import NGSQuality, NGSQualityEncoding, QualityMatrix, NGSQualityBaseOffset
from seqtools.command import Command

import qualities
import fastq
import fasta
import qseq
import bowtie

# __all__ = ['fasta', 'solid']

