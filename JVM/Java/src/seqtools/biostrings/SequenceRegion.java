package seqtools.biostrings;

class SequenceRegion<T> implements Sequence {
  
  protected T sequence;
  protected long start;
  protected long end;
  
  protected SequenceRegion() {}
  public SequenceRegion(T seq, long start, long end) {
    if (end < 0 || start < 0) {
      throw new IllegalArgumentException("start and end must be postive");
    }
    if (end < start) {
      throw new IllegalArgumentException("end must be greater than start");
    }
    if (end >= seq.length()) {
      throw new IllegalArgumentException("end position longer than sequence");
    }
    this.sequence = seq;
    this.start = start;
    this.end = end;
  }
  
  @Override
  public long length() {
    return this.end - this.start + 1;
  }
  
  @Override
  public Iterable<T> getSequence() {
    return this.sequence;
  }
}