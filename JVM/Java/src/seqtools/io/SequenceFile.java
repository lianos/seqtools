package seqtools.io;

import java.io.File;

public interface SequenceFile {
  public boolean hasNext();
  public Sequence next() throws SequenceFormatException;
  public boolean isColorspace();
  public String name();
  public int getPercentComplete();
  public File getFile();
}