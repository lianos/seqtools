package seqtools.biostrings;

public abstract class Sequence {
  
  protected byte[] sequence;
  protected String name;
  
  public Sequence(String name, byte[] seq) {
    this.name = name;
    this.sequence = seq;
  }
  
  public view(long start, long end) {
    return new SequenceRegion(this, start, end);
  }
  
  public long length() {
    this.sequence.length;
  }
}
