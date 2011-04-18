package seqtools.io;

/**
 * Design taken from FastQC BAMFile
 */
import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class BamFile implements AlignmentFile {
    
    private File file;
    private boolean onlyMapped;
    private long fileSize = 0;
    private long recordSize = 0;
    private int count = 0;
    private int rawCount = 0;
    private SAMFileReader br;
    private String name;
    private Sequence nextSequence = null;
    Iterator<SAMRecord> it;
    
    
    protected BAMFile (File file, boolean onlyMapped) throws
            FileFormatException, IOException {
        this.file = file;
        fileSize = file.length();
        name = file.getName();
        this.onlyMapped = onlyMapped;
        SAMFileReader.setDefaultValidationStringency(
            SAMFileReader.ValidationStringency.SILENT
        );
        
        br = new SAMFileReader(file);
        it = br.iterator();
        readNext();
    }
    
    
    private void readNext() {
        SAMRecord record;
        while (true) {
            if (!it.hasNext()) {
                nextSequence = null;
                return;
            }
            record = it.next();
            ++rawCount;
            // We skip over entries with no mapping if that's what the user 
            // asked for
            if (onlyMapped && record.getReadUnmappedFlag()) {
                continue;
            }
            else {
                break;
            }
        }
        ++count;
        if (recordSize == 0) {
            recordSize = (record.getReadLength()*2)+150;
            if (br.isBinary()) {
                recordSize /= 4;
            }
        }
        String sequence = record.getReadString();
        String qualities = record.getBaseQualityString();
        // BAM/SAM files always show sequence relative to the top strand of
        // the mapped reference so if this sequence maps to the reverse strand
        // we need to reverse complement the sequence and reverse the qualities
        // to get the original orientation of the read.
        if (record.getReadNegativeStrandFlag()) {
            sequence = reverseComplement(sequence);
            qualities = reverse(qualities);
        }
        nextSequence = new Sequence(this, sequence, qualities, record.getReadName());
    } // readNext
    
}