package seqtools.aligner.bowtie

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Runs through an bowtie-generated SAM file to filter out duplicate reads.
 * 
 * For now we it keeps only reads that have less than <code>N</code> alignments
 * in the top tier.
 * 
 * Multiple alignments are reported <emph>sequentially</emph> and in order
 * of quality.
 */
public class UniqueAlignmentFilter {
    
    AlignmentFile file;
}
